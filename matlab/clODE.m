classdef clODE < cppclass & matlab.mixin.SetGet
    % clODE(prob, stepper=rk4, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
    
    %TODO: Using set/get mixin causes setters to be called during
    %constructor; unnecessary duplicate calls to C++ functions. 
    
    properties
        
        stepper='dopri5'
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
                stepper='dopri5';
            end
            
            devices=queryOpenCL(); %throws error if no opencl 
            %device selection: parse inputs
            if  ~exist('selectedDevice','var')
                selectedDevice=[];
            elseif selectedDevice>length(devices)
                error('Device index specifed is greater than the number of OpenCL devices present')
            end
            %auto-selection of devices: gpu>cpu
            if isempty(selectedDevice)
                selectedDevice=clODE.autoselectDevice(devices);
            end
            
            args{1}=prob;
            args{2}=stepper;
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
            obj.sp=clODE.defaultSolverParams(); %default solver params
        end
        
        % new and delete are inherited
        
        %set a new problem - must initialize again!
        function setNewProblem(obj, prob)
            if ~strcmp(prob,obj.prob)
                obj.prob=prob;
                obj.cppmethod('setnewproblem', prob);
            end
        end
        
        %set a new time step method - must initialize again!
        function set.stepper(obj, newStepper)
            if ~strcmp(newStepper,obj.stepper)
                obj.stepper=newStepper;
                obj.cppmethod('setstepper', newStepper);
            end
        end
        
        %set single precision true/false - must initialize again! 
        %TODO: Need to run ode2cl again!!! Alt: just generate both and swap filenames...
        function set.precision(obj, newPrecision)
            if ~strcmp(newPrecision,obj.precision)
                switch lower(newPrecision)
                    case {'single'}
                        obj.precision='single';
                        clSinglePrecision=true;
                    case {'double'}
                        obj.precision='double';
                        clSinglePrecision=false;
                    otherwise
                        error('Precision must be set to ''single'' or ''double''')
                end
%                 [~,obj.prob]=ode2cl(obj.prob.file,[],clSinglePrecision);
                obj.cppmethod('setprecision', clSinglePrecision);
            end
        end
        
        %set a new OpenCL context - must initialize again!
        function setOpenCL(obj, newDevice)
            if newDevice~=obj.selectedDevice && newDevice<=length(obj.devices)
                obj.selectedDevice=newDevice;
                obj.cppmethod('setopencl', obj.devices(newDevice).platformID, obj.devices(newDevice).deviceID);
            end
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
        function transient(obj, tspan)
            if exist('tspan','var')
                obj.settspan(tspan);
            end
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
        
        function stepperNames=getAvailableSteppers(obj)
            stepperNames=obj.cppmethod('getsteppernames');
        end
        function programString=getProgramString(obj)
            prog=obj.cppmethod('getprogramstring');
            programString=sprintf('%s',prog{1});
        end
        function printStatus(obj)
            obj.cppmethod('printstatus');
        end
    end
    
    
    %static helper methods
    methods (Static=true)
        
        function sp=defaultSolverParams()
            sp.dt=.1;
            sp.dtmax=100.00;
            sp.abstol=1e-6;
            sp.reltol=1e-3;
            sp.max_steps=1000000;
            sp.max_store=10000; %allocated number of timepoints: min( (tf-t0)/(dt*nout)+1 , sp.max_store)
            sp.nout=1;
        end
        
        function selectedDevice=autoselectDevice(devices)
            selectedDevice=[];
%             if  isempty(selectedDevice)
%                 selectedDevice=find({devices(:).type}=="accel");
%             end
            if  isempty(selectedDevice)
                selectedDevice=find({devices(:).type}=="GPU");
            end
            if  isempty(selectedDevice)
                selectedDevice=find({devices(:).type}=="CPU");
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

