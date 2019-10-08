classdef clODE < cppclass & matlab.mixin.SetGet
    % clODE(prob, stepper=rk4, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
    
    %TODO: this should have a method for default solverparams
    
    properties
        
        stepper='rk4'
        precision='single'
        
%         cl_vendor='any'
%         cl_deviceType='default'
        devices %array of structs describing OpenCL compatible devices
        selectedDevice=1
        
        prob
        
        nPts
        P
        X0
        Xf
        
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
        
        function obj = clODE(arg1, precision, selectedDevice, stepper, mexFilename, extraArgs)
            
            if nargin==0
                error('Problem info struct is a required argument')
            end
            
            if  ~exist('precision','var')||isempty(precision)
                precision='single';
            end
            clSinglePrecision=true;
            if precision=="double", clSinglePrecision=false;end
            
            if ischar(arg1)
                [~,prob]=ode2cl(arg1,[],clSinglePrecision);
            elseif isstruct(arg1)
                prob=arg1;
            else
%                 if nargin==0 || isempty(source)
%                     [name,path]=uigetfile('.ode','Select an ODE file');
%                     if ~ischar(name)
%                         disp('File selection canceled, quitting...')
%                         return
%                     end
%                     source=fullfile(path,name);
%                 end
            end
            
            if  ~exist('stepper','var')||isempty(stepper)
                stepper='dorpri5';
            end
            
            devices=queryOpenCL(); %default: first device
            if  ~exist('selectedDevice','var')||isempty(selectedDevice)
                selectedDevice=1;
            end
            
            args{1}=prob;
            args{2}=clODE.getStepperEnum(stepper);
            args{3}=clSinglePrecision;
            args{4}=devices(selectedDevice).platformID;
            args{5}=devices(selectedDevice).deviceID;

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
            obj.precision=precision;
            obj.devices=devices;
            obj.selectedDevice=selectedDevice;
        end
        
        % new and delete are inherited
        
        %set a new problem - must initialize again!
        function set.prob(obj, prob)
            obj.prob=prob;
            obj.cppmethod('setnewproblem', prob);
        end
        
        %set a new time step method - must initialize again!
        function set.stepper(obj, newStepper)
            stepperEnum=clODE.getStepperEnum(newStepper);
            obj.stepper=newStepper;
            obj.cppmethod('setstepper', stepperEnum);
        end
        
        %set single precision true/false - must initialize again! 
        %TODO::: Need to run ode2cl again???
        function set.precision(obj, newPrecision)
            clSinglePrecision=true; %default to single
            obj.precision='single';
            switch lower(newPrecision)
                case {'single',1}
                case {'double',2}
                    obj.precision='double';
                    clSinglePrecision=false;
                otherwise
                    warning('Precision must be set to ''single'' or ''double''. Using single precision.')
            end
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
        
        function shiftTspan(obj)
            obj.cppmethod('shifttspan');
        end
        
        function shiftX0(obj)
            obj.cppmethod('shiftx0');
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
        
        function Xf=getXf(obj)
            Xf=obj.cppmethod('getxf');
            Xf=reshape(Xf,obj.nPts,obj.prob.nVar);
            obj.Xf=Xf;
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
                case {'dorpri5','dorpri'}
                    stepperInt=5;
                otherwise
                    warning('Unrecognized stepper method name. Using default: dorpri5')
                    stepperInt=5;
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

