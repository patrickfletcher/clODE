classdef clODEtrajectory<clODE & matlab.mixin.SetGet
    % clODEtrajectory(prob, stepper=rk4, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
    
    %TODO: return only valid time points for each trajectory - cell array??
    
    properties
        nStored
        
        t
        x
        dx
        aux
        
    %inherits from clODE:
%     prob
%     stepper
%     clSinglePrecision
%     cl_vendor
%     cl_deviceType
%     
%     nPts
%     P
%     X0
%     Xf
%     auxf
%     sp
%     tspan
    end
    
    methods
       
        function obj = clODEtrajectory(arg1, precision, selectedDevice, stepper, mexFilename)
            %hack to get correct mexfile for classes derived from this one.
            %I want the subclass to get all the methods contained here, but
            %it needs to use a mex function that unfortunately has to
            %repeat base class method dispatch code.
            
            if  ~exist('precision','var')||isempty(precision)
                precision=[]; %default handled in clODE.m
            end
            
            if  ~exist('selectedDevice','var')||isempty(selectedDevice)
                selectedDevice=[]; %default handled in clODE.m
            end
            
            if  ~exist('stepper','var')||isempty(stepper)
                stepper=[]; %default handled in clODE.m
            end
            
            if ~exist('mexFilename','var')
                mexFilename='clODEtrajectorymex';
            end
            obj@clODE(arg1, precision, selectedDevice, stepper, mexFilename);
        end
        
        
        %overloads to fetch data if desired
        function trajectory(obj)
            obj.cppmethod('trajectory');
            if any(obj.nStored==obj.sp.max_store)
                warning('Maximum storage reached')
            end
        end
        
        function t=getT(obj)
            t=obj.cppmethod('gett'); t=t(:); %force column
            t=reshape(t,obj.nPts,obj.sp.max_store)';
            obj.t=t;
        end
        
        function x=getX(obj)
            x=obj.cppmethod('getx');
            x=reshape(x,obj.nPts,obj.prob.nVar,obj.sp.max_store);
            x=permute(x,[3,2,1]);
            obj.x=x;
        end
        
        function dx=getDx(obj)
            dx=obj.cppmethod('getdx');
            dx=reshape(dx,obj.nPts,obj.prob.nVar,obj.sp.max_store);
            dx=permute(dx,[3,2,1]);
            obj.dx=dx;
        end
        
        function aux=getAux(obj)
            nAux=max(obj.prob.nAux,1); %make a dummy aux trajectory if nAux=0
            aux=nan(obj.sp.max_store,nAux,obj.nPts);
            obj.aux=aux;
            if obj.prob.nAux>0
            aux=obj.cppmethod('getaux');
            aux=reshape(aux,obj.nPts,obj.prob.nAux,obj.sp.max_store);
            aux=permute(aux,[3,2,1]);
            obj.aux=aux;
            end
        end

        function nStored=getNstored(obj)
            nStored=obj.cppmethod('getnstored');
            obj.nStored=nStored+1;
        end
        
    end
    
end