classdef clODEtrajectory<clODE
    % clODEtrajectory(prob, stepper=rk4, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
    
    %TODO: return only valid time points for each trajectory - cell array??
    
    properties
        nSteps
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
%     auxf
%     sp
%     tspan
    end
    
    methods
       
        function obj = clODEtrajectory(prob, stepper, clSinglePrecision, cl_vendor, cl_deviceType, mexFilename)
            %hack to get correct mexfile for classes derived from this one.
            %I want the subclass to get all the methods contained here, but
            %it needs to use a mex function that unfortunately has to
            %repeat base class method dispatch code.
            if ~exist('mexFilename','var')
                mexFilename='clODEtrajectorymex';
            end
            obj@clODE(prob, stepper, clSinglePrecision, cl_vendor, cl_deviceType, mexFilename);
        end
        
        
        %overloads to fetch data if desired
        function trajectory(obj)
            obj.cppmethod('trajectory');
            obj.getNsteps();
        end
        
            
        function nSteps=getNsteps(obj)
            nSteps=obj.cppmethod('getsteps');
            obj.nSteps=nSteps;
        end
        function t=getT(obj)
            t=obj.cppmethod('gett'); t=t(:); %force column
            t=reshape(t,obj.nPts,obj.nSteps)';
            obj.t=t;
        end
        
        function x=getX(obj)
            x=obj.cppmethod('getx');
            x=reshape(x,obj.nPts,obj.prob.nVar,obj.nSteps);
            x=permute(x,[3,2,1]);
            x=squeeze(x);
            obj.x=x;
        end
        
        function dx=getDx(obj)
            dx=obj.cppmethod('getdx');
            dx=reshape(dx,obj.nPts,obj.prob.nVar,obj.nSteps);
            dx=permute(dx,[3,2,1]);
            dx=squeeze(dx);
            obj.dx=dx;
        end
        
        function aux=getAux(obj)
            aux=obj.cppmethod('getaux');
            aux=reshape(aux,obj.nPts,obj.prob.nAux,obj.nSteps);
            aux=permute(aux,[3,2,1]);
            aux=squeeze(aux);
            obj.aux=aux;
        end

        function nStored=getNstored(obj)
            nStored=obj.cppmethod('getstored');
            obj.nStored=nStored;
        end
        
    end
    
end