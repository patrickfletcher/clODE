function [pointset,num_points,ranged_dim,grid_dim]=generate_pointset(type, varargin)
% Collection of methods for generating point sets. The first argument
% defines the type of point set to create. Each type requires additional
% input arguments, detailed below. The result is matrix containing a
% point's coordinates in each row.
%
% USAGE
%   X=generate_pointset('grid',lb,ub,grid_dim)
%     lb, ub: lower and upper bounds for the box containing points. If
%     lb=ub, that value is used for all points.
%     grid_dim: scalar=symmetric grid, vector=grid with given number of
%     points in each dimension.
%
%   X=generate_pointset('random',lb,ub,num_points) - uses "rand"
%
%   X=generate_pointset('latin',lb,ub,num_points) - latin hypercube sampling
%
%   X=generate_pointset('sobol',lb,ub,num_points) - Sobol quasi-random
%   sequence
%
%   X=generate_pointset('halton',lb,ub,num_points) - Halton quasi-random
%   sequence
%
%   X=generate_pointset('gaussian',origin_point,sigma,num_points,ranged_dims)
%      -uses randn with mean at origin_point and covariance matrix sigma.
%      Optionally supply ranged_dims (list of which dimensions to use) to
%      generate points in a subset of coordinates
%
%   X=generate_pointset('repeated',origin_point,num_points) - simple repeat

%   X=generate_pointset('coordinate_planes',origin_point,lb,ub,type,typearg)
%     Recursively calls GeneratePointSet for each combination of two
%     coordinates, generating for each a plane of 

% TODO: add support for ('name','value') pairs as more arguments to
% control the specifics of each type of point set. e.g. rng seed, latin
% hypercube refinement, quasi-random pointset controls (skip, leap,
% scramble)

box_pointset=false;
grid_dim=[];
ranged_dim=[];
num_points=[];
pointset=[];

%check and unpack inputs
switch type
    case 'grid'
        %extra inputs must include: lb, ub, grid_dimensions
        lb=varargin{1};
        ub=varargin{2};
        grid_dim=varargin{3};
        
        if length(ub)~=length(lb)
            error('lb and ub must have the same dimension');
        end
        
        dim=length(lb);
        ranged_dim=(ub-lb)>0;
        num_ranged=nnz(ranged_dim);
        n_grid=numel(grid_dim);
        
        if n_grid==1
            grid_dim=repmat(grid_dim,1,num_ranged);
        elseif n_grid>1 && n_grid~=num_ranged
            error('The specified grid dimensions don''t match the number of variables with (ub-lb)>0');
        end
        
        box_pointset=true;
        
    case {'random','latin','halton','sobol'}
        %extra inputs must include: lb, ub, num_points
        lb=varargin{1};
        ub=varargin{2};
        num_points=varargin{3};
        
        dim=length(lb);
        ranged_dim=(ub-lb)>0;
        num_ranged=nnz(ranged_dim);
        
        if length(ub)~=dim
            error('lb and ub must have the same dimension');
        end
        
        box_pointset=true;
        
    case 'gaussian'
        %extra inputs must include: origin_point, num_points, ranged_dim
        origin_point=varargin{1};
        sig=varargin{2}; %for now assuming constant...
        num_points=varargin{3};
        
        if length(varargin)<4
            ranged_dim=[];
        else
            ranged_dim=varargin{4};
        end
        
        dim=length(origin_point);
        if isempty(ranged_dim)
            ranged_dim=ones(1,dim);
        else
            
        ranged_dim=varargin{4};
        end
        num_ranged=nnz(ranged_dim);
        
    case 'repeated'
        origin_point=varargin{1};
        num_points=varargin{2};
        
        dim=length(origin_point);
        ranged_dim=ones(1,dim);
        num_ranged=nnz(ranged_dim);
end

%build the points
switch type

    case 'grid'
        
        num_points=prod(grid_dim);
        X=makeNDgrid(grid_dim);
        
        
        %NOTE: This is equivalent! can use reshape to go back to grid
        
%         per=[6,10,20,30,60,90,120,240];
%         x=linspace(0,10,16);
%         y=linspace(1,5,16);
%         [X,Y,P]=ndgrid(x,y,per);
%         XX=[X(:),Y(:),P(:)];
        

    case 'random'
        
        X=rand(num_points,num_ranged);

    case 'latin'
        X=lhsdesign(num_points,num_ranged,'criterion','none');
%         X=lhsdesign(num_points,num_ranged,'criterion','maximin','iterations',5); %Maximize minimum distance between points
%         X=lhsdesign(num_points,num_ranged,'criterion','correlation','iterations',5); %Reduce correlation

    case 'halton'
        %RR2 scramble on Halton is a deterministic algorithm
        p=haltonset(num_ranged,'Skip',1e3,'Leap',1e2);
        p=scramble(p, 'RR2');
        X=net(p,num_points);
        
    case 'sobol'
        p=sobolset(num_ranged,'Skip',1e3,'Leap',1e2);
        p=scramble(p, 'MatousekAffineOwen');
        X=net(p,num_points);
        
    case 'gaussian'
        mu=origin_point(ranged_dim);
        
        %from the matlab help: using covariance matrix?
%         sig = [1 0.5; 0.5 2];
%         R = chol(sigma);
%         z = repmat(mu,10,1) + randn(10,2)*R

        X=repmat(mu,num_points,1)+ randn(num_points,num_ranged)*sig;
        
    case 'repeated'
        %for convenience, repmat of origin_points.
        X=repmat(origin_point,num_points,1);
        
    otherwise
        error('unknown point type requested');
end

%put the ranged values into the correct columns of the output array
pointset=ones(num_points,dim);
pointset(:,ranged_dim)=X;
        
%some of the types were generated in the unit cube, must be rescaled:
if box_pointset
    for i=1:dim
        pointset(:,i)=lb(i)+(ub(i)-lb(i))*pointset(:,i);
    end
end

end


function X=makeNDgrid(grid_size)
    %make a grid of prod(nGrid) pts with spacing in each dimension given in
    %each element of nGrid
    
    grid_dim=nnz(grid_size>1);
    nPts=prod(grid_size);
    X=ones(nPts,grid_dim);
    
    dx=1./(grid_size-1); %unit interval in each dimension, will rescale after
    
    %use the linear index idea from ind2sub
    cum=[1,cumprod(grid_size(1:end-1))];
    for i=1:nPts
        idx=i;
        vi=0;
        vj=0;
        for k=grid_dim:-1:1
            vi=rem(idx-1,cum(k))+1;
            vj=(idx-vi)/cum(k);
            idx=vi;
            
            X(i,k)=vj*dx(k);
        end
    end
    
end
