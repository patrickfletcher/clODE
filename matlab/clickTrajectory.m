function clickTrajectory(src,evnt,clo,Grid,click_pars,markerColor,X0grid,Pgrid)
%A mouse click inside the axes will trigger a parameter
%selection crosshair. A second click selects a parameter
%combination to simulate, and resulting simulation is displayed
%in a new figure window.

%could find closest solved instance, get x0 for beginning of
%feature detection interval


if ~exist('click_pars','var') || isempty(click_pars)
    click_pars.vars=clo.prob.varNames(1); 
    click_pars.nClick=1;
    click_pars.fig=[];
    click_pars.same_ylims=0;
end
nClick=click_pars.nClick;
vars=click_pars.vars;
trajFigID=click_pars.fig;

if ischar(vars), vars={vars}; end
if nClick>6, error('I refuse to do more than 6 clicks'); end

if ~exist('markerColor','var'), markerColor='k';end
doX0search=false;
if exist('X0grid','var')&&exist('Pgrid','var'), doX0search=true;end

tspan=clo.tspan;
x0=repmat(clo.prob.x0(:)',nClick,1);
p=repmat(clo.prob.p0(:)',nClick,1);

% determine whether desired variables are var/aux %TODO: homogenize by
% combining var+aux into one list
varIsAux=false(size(vars));
for v=1:length(vars)
    tmpIx=find(clo.prob.varNames==string(vars{v}));
    if isempty(tmpIx)
        tmpIx=find(clo.prob.auxNames==string(vars{v}));
        if isempty(tmpIx)
            error(['unknown variable: ' vars{v}]);
        else
            varIsAux(v)=true;
            varIx(v)=tmpIx;
        end
    else
        varIx(v)=tmpIx;
    end
end

markers='osdp^v';

oldMarks=findobj('tag','marks');
if ~isempty(oldMarks)
    delete(oldMarks);
end

for n=1:nClick
    [px,py]=ginput(1);
    p(n,Grid.xix)=px;
    p(n,Grid.yix)=py;
    h(n)=line(px,py,'marker',markers(n),'linestyle','none','color',markerColor,'linewidth',1,'tag','marks');
    
    if doX0search
        pix = knnsearch(Pgrid,[px,py]);
        x0(n,:)=X0grid(pix,:);
    end
end


clo.setProblemData(x0,p);
clo.settspan(tspan);
% clo.transient();

if ~exist('trajFigID','var')||isempty(trajFigID)
    trajFigID=figure();
end
hf=figure(trajFigID);
clf reset

hf.KeyPressFcn=@keypress;

clf
[X,T,AUX]=integrate();

%TODO: if this is a chart class, would setup axes/line objects, then just
%update data+lims
gap=[0.05,0.03];marg_h=[0.1,0.05]; marg_w=[0.075,0.075];
ax=tight_subplot(nClick,1,[],gap,marg_h,marg_w);

plotTrajectories(X,T,AUX);

    function [X,T,AUX]=integrate()
        clo.trajectory();
        clo.getXf();
        xx=clo.getX();
        aux=clo.getAux();
        tt=clo.getT()*clo.tscale;
        nStored=clo.getNstored();
        
        %support adaptive timestepping: each trajectory may have different
        %number of steps. array of structs seems conceptually better
%         traj=struct('t',[],'x',[],'aux',[]);
        for i=1:nClick
            T{i}=tt(1:nStored(i),i);
            X{i}=xx(1:nStored(i),:,i);
            AUX{i}=aux(1:nStored(i),:,i);
        end
    end

    function plotTrajectories(X,T,AUX)
        
        %extract data to plot
        lymin=inf; lymax=-inf;
        rymin=inf; rymax=-inf;
        traj=struct('t',[],'x',[],'xname',{},'tlim',[],'xlo',[],'xhi',[]);
        for i=1:nClick
            traj(i).t=T{i};
            for vv=1:length(varIx)
                if varIsAux(vv)
                    traj(i).x(:,vv)=AUX{i}(:,varIx(vv));
                    traj(i).xname(vv)=clo.prob.auxNames(varIx(vv));
                else
                    traj(i).x(:,vv)=X{i}(:,varIx(vv));
                    traj(i).xname(vv)=clo.prob.varNames(varIx(vv));
                end
            end
            %individual data limits
            traj(i).tlim=[traj(i).t(1), traj(i).t(end)];
            traj(i).xlo=min(traj(i).x,[],1);
            traj(i).xhi=max(traj(i).x,[],1);
        end
        %global data limits
        XLO=min(cat(1,traj(:).xlo),[],1);
        XHI=max(cat(1,traj(:).xhi),[],1);
        TLIMS=cat(1, traj(:).tlim);
        TLIM=[min(TLIMS(:,1)), max(TLIMS(:,2))];
        
        for i=1:nClick
            
            axes(ax(i)); cla reset
            
            %left axis
            plot(traj(i).t,traj(i).x(:,1));
            ylabel(traj(i).xname(1));
            
            xlim(TLIM);
            if click_pars.same_ylims
                ylim([XLO(1),XHI(1)]);
            end
%             ylim([min(YLIM(1),xmin-0.01*abs(xmin)), max(YLIM(2),xmax+0.01*abs(xmax))]);

            if numel(varIx)==2
                yyaxis right
                plot(traj(i).t,traj(i).x(:,2));
                xlim(TLIM);
                if click_pars.same_ylims
                    ylim([XLO(2),XHI(2)]);
                end
%                 ylim([min(YLIM(1),xmin-0.01*abs(xmin)), max(YLIM(2),xmax+0.01*abs(xmax))]);
                ylabel(traj(i).xname(2));
            end
            
            box off
            
            if i<nClick
                ax(i).XTickLabel=[];
            else
                ax(i).XTickLabelMode='auto';
                xlabel(['t' ' (' clo.tunits ')']);
            end
            
            XLIM=xlim();
            YLIM=ylim();
            line((XLIM(1)+XLIM(end))*0.02,YLIM(2)*0.9,'marker',markers(i),...
                'linestyle','none','color',markerColor,'linewidth',1,'tag','marks');
            title([clo.prob.parNames{Grid.xix} '=' num2str(p(i,Grid.xix), '%.4g') ', '...
                clo.prob.parNames{Grid.yix} '=' num2str(p(i,Grid.yix),'%.4g')]);
            
        end
        
        linkaxes(ax,'x')
    end

    %interactivity
    function keypress(src,evt)
        switch(evt.Key)
            case 'c'
                clo.shiftTspan();
                clo.shiftX0();
                [xx,tt,aux]=integrate();
                for i=1:nClick
                    T{i}=[T{i};tt{i}];
                    X{i}=[X{i};xx{i}];
                    AUX{i}=[AUX{i};aux{i}];
                end
                plotTrajectories(X,T,AUX);
                
            case 'g' %go
                clo.settspan(tspan);
                [X,T,AUX]=integrate();
                plotTrajectories(X,T,AUX);
                
            case 'l' %last
                clo.settspan(tspan);
                clo.shiftX0();
                [X,T,AUX]=integrate();
                plotTrajectories(X,T,AUX);
                
            case 'r' %randomize ICs
                x0lb=[clo.prob.var.lb];
                x0ub=[clo.prob.var.ub];
                x0=x0lb+rand(nClick,length(x0lb)).*(x0ub-x0lb);
                clo.settspan(tspan);
                clo.setX0(x0(:));
                [X,T,AUX]=integrate();
                plotTrajectories(X,T,AUX);
                
            case 's' %shift
                clo.shiftTspan();
                clo.shiftX0();
                [X,T,AUX]=integrate();
                plotTrajectories(X,T,AUX);
        end
    end
end