classdef clODE < cppclass & matlab.mixin.SetGet
    % clODE(prob, stepper=rk4, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
    
    %TODO: this should have a method for default solverparams
    
    properties
        
        stepper='rk4'
        clSinglePrecision=true
        cl_vendor='any'
        cl_deviceType='default'
        
        prob
        
        nPts
        P
        X0
        
        sp
        tspan
        
        tscale=1 %should go in IVP  
        tunits='';
    end
    
    properties (Dependent)
        %nPts
    end
    
    properties (Access = private)
    end
    
    
    %C++ class interaction
    methods
        
        function obj = clODE(prob, stepper, clSinglePrecision, cl_vendor, cl_deviceType, mexFilename, extraArgs)
            
            if nargin==0
                error('Problem info struct is a required argument')
            end
            
            args{1}=prob;
            
            if  ~exist('stepper','var')||isempty(stepper)
                stepper='rungekutta4';
            end
            args{2}=clODE.getStepperEnum(stepper);
            
            if  ~exist('clSinglePrecision','var')||isempty(clSinglePrecision)
                clSinglePrecision=true;
            end
            args{3}=clSinglePrecision;
            
            if  ~exist('cl_vendor','var')||isempty(cl_vendor)
                cl_vendor='any';
            end
            args{4}=clODE.getVendorEnum(cl_vendor);
            
            if  ~exist('cl_deviceType','var')||isempty(cl_deviceType)
                cl_deviceType='default';
            end
            args{5}=clODE.getDeviceTypeEnum(cl_deviceType);
            
            %hack to get correct mexfile for classes derived from this one.
            %I want the subclass to get all the methods contained here, but
            %it needs to use a mex function that unfortunately has to
            %repeat base class method dispatch code.
            if ~exist('mexFilename','var')||isempty(mexFilename)
                mexFilename='clODEmex';
            end
            
            if exist('extraArgs','var')
                args=[args,extraArgs];
            end
            
            obj@cppclass(mexFilename,args{:});
            
            obj.prob=prob;
            obj.stepper=stepper;
            obj.clSinglePrecision=clSinglePrecision;
            obj.cl_vendor=cl_vendor;
            obj.cl_deviceType=cl_deviceType;
        end
        
        % new and delete are inherited
        
        %set a new problem - must initialize again!
        function set.prob(obj, prob)
            obj.prob=prob;
            obj.cppmethod('setnewproblem', prob);
        end
        
        %set a new time step method - must initialize again!
        function set.stepper(obj, newStepper)
            obj.stepper=newStepper;
            obj.cppmethod('setstepper', clODE.getStepperEnum(newStepper));
        end
        
        %set single precision true/false - must initialize again!
        function setPrecision(obj, clSinglePrecision)
            obj.clSinglePrecision=clSinglePrecision;
            obj.cppmethod('setprecision', clSinglePrecision);
        end
        
        %set a new OpenCL context - must initialize again!
        function setOpenCL(obj, newVendor, newDeviceType)
            obj.cl_vendor=newVendor;
            obj.cl_deviceType=newDeviceType;
            vendorInt=clODE.getVendorEnum(newVendor);
            deviceTypeInt=clODE.getDeviceTypeEnum(newDeviceType);
            obj.cppmethod('setopencl', vendorInt, deviceTypeInt);
        end
        
        %initialize builds the program and sets data needed to run
        %simulation in one call
        function initialize(obj, tspan, X0, P, sp)
            obj.cppmethod('initialize', tspan, X0(:), P(:), sp);
            obj.tspan=tspan;
            obj.X0=X0;
            obj.P=P;
            obj.sp=sp;
            obj.nPts=numel(X0)/obj.prob.nVar;
        end
        
        %Set X0 and P together if trying to change nPts
        function setProblemData(obj, X0, P)
            obj.X0=X0;
            obj.P=P;
            obj.nPts=numel(X0)/obj.prob.nVar;
            obj.cppmethod('setproblemdata', X0(:), P(:));
        end
        
        
        function seedRNG(obj, mySeed)
            if ~exist('mySeed','var')
                obj.cppmethod('seedrng');
            else
                obj.cppmethod('seedrng', mySeed);
            end
        end
        
        
        function settspan(obj, tspan)
            obj.tspan=tspan;
            obj.cppmethod('settspan', tspan);
        end
        
        %nPts cannot change here
        function setX0(obj, X0)
            testnPts=numel(X0)/obj.prob.nVar;
            if testnPts==obj.nPts
                obj.X0=X0;
                obj.cppmethod('setx0', X0(:));
            else
                error('Size of X0 is incorrect');
            end
        end
        
        %nPts cannot change here
        function setP(obj, P)
            testnPts=numel(P)/obj.prob.nPar;
            if testnPts==obj.nPts
                obj.P=P;
                obj.cppmethod('setpars', P(:));
            else
                error('Size of P is incorrect');
            end
        end
        
        function setsp(obj, sp)
            obj.sp=sp;
            obj.cppmethod('setsolverpars', sp);
        end
        
        %matlab object becomes de-synchronized from GPU arrays when calling
        %simulation routines. User must trigger data fetch from GPU
        function transient(obj)
            obj.cppmethod('transient');
        end
        
        function tspan=getTspan(obj)
            tspan=obj.cppmethod('gettspan');
            obj.tspan=tspan;
        end
        
        function X0=getX0(obj)
            X0=obj.cppmethod('getx0');
            X0=reshape(X0,obj.nPts,obj.prob.nVar);
            obj.X0=X0;
        end
        
        function auxf=getAuxf(obj)
            auxf=obj.cppmethod('getauxf');
            auxf=reshape(auxf,obj.nPts,obj.prob.nAux);
        end
        
    end
    
    
    %static helper methods
    methods (Static=true)
        
        function stepperInt=getStepperEnum(steppername)
            switch lower(steppername)
                case {'euler'}
                    stepperInt=0;
                case {'heun','modeuler'}
                    stepperInt=1;
                case {'rk4','runge','rungekutta4'}
                    stepperInt=2;
                case {'heuneuler'}
                    stepperInt=3;
                case {'bs23','bogackishampine'}
                    stepperInt=4;
                case {'dorpri5'}
                    stepperInt=5;
                otherwise
                    warning('Unrecognized stepper method name. Using default: RK4')
                    stepperInt=2;
            end
        end
        
        function vendorInt=getVendorEnum(cl_vendor)
            switch lower(cl_vendor)
                case {'any'}
                    vendorInt=0;
                case {'nvidia'}
                    vendorInt=1;
                case {'amd'}
                    vendorInt=2;
                case {'intel'}
                    vendorInt=3;
                otherwise
                    warning('Unrecognized vendor name. Using default: any vendor')
                    vendorInt=0;
            end
        end
        
        function deviceTypeInt=getDeviceTypeEnum(cl_deviceType)
            switch lower(cl_deviceType)
                case {'default'}
                    deviceTypeInt=0;
                case {'cpu'}
                    deviceTypeInt=1;
                case {'gpu'}
                    deviceTypeInt=2;
                case {'accel','accelerator'}
                    deviceTypeInt=3;
                case {'all'}
                    deviceTypeInt=4;
                otherwise
                    warning('Unrecognized device type name. Using default device type')
                    deviceTypeInt=0;
            end
        end
        
    end
    
end

