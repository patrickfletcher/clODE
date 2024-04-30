classdef clODE < cppclass
    % clODE(prob, stepper=rk4, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
    
    %TODO: Using set/get mixin causes setters to be called during
    %constructor; unnecessary duplicate calls to C++ functions. 
    
    properties
        
%         cl_vendor='any'
%         cl_deviceType='default'
        devices %array of structs describing OpenCL compatible devices
        selectedDevice=[]
        
        stepper='rk4'
        precision='single'
        
        prob
        
        nPts
        P
        X0
        Xf
        
        sp
        tspan
        
        clBuilt=false
        clInitialized=false
        
        tscale=1 %should go in IVP  
        tunits=''
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
            
            if ischar(arg1)||isstring(arg1)
                arg1=char(arg1);
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
            if ~exist('selectedDevice','var')
                selectedDevice=[];
            elseif selectedDevice>length(devices)
                error('Device index specifed is greater than the number of OpenCL devices present')
            end
            %auto-selection of devices: gpu>cpu
            if isempty(selectedDevice)
                selectedDevice=clODE.autoselectDevice(devices);
            end

            % purify the problem struct
            newprob = clODE.purifyProblemStruct(prob);
            
            args{1}=newprob;
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
            
            %TODO: mex call constructing the object is slow!!! why?
            obj@cppclass(mexFilename,args{:});
            
            obj.prob=newprob;
            obj.devices=devices;
            obj.selectedDevice=selectedDevice;
            
            obj.stepper=stepper;
            obj.precision=precision;
            obj.sp=clODE.defaultSolverParams(); %default solver params (device transfer during init)
%             obj.buildCL();
        end
        
        % new and delete are inherited
        
        %set a new problem - must initialize again!
        function setProblemInfo(obj, newprob)
            newprob = clODE.purifyProblemStruct(newprob);
            if ~strcmp(newprob,obj.prob)
                obj.cppmethod('setProblemInfo', newprob);
                obj.prob=newprob;
                obj.clBuilt=false;
                obj.clInitialized=false;
            end
        end
        
        %set a new time step method - must initialize again!
        function setStepper(obj, newStepper)
            if ~strcmp(newStepper,obj.stepper)
                obj.stepper=newStepper;
                obj.cppmethod('setstepper', newStepper);
                obj.clBuilt=false;
                obj.clInitialized=false;
            end
        end
        
        %set single precision true/false - must initialize again! 
        %ode2cl generates a file with "realtype"
        function setPrecision(obj, newPrecision)
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
                obj.cppmethod('setprecision', clSinglePrecision);
                obj.clBuilt=false;
                obj.clInitialized=false;
            end
        end
        
        %set a new OpenCL context - must initialize again!
        function setOpenCL(obj, newDevice)
            if newDevice~=obj.selectedDevice && newDevice<=length(obj.devices)
                obj.selectedDevice=newDevice;
                obj.cppmethod('setopencl', obj.devices(newDevice).platformID, obj.devices(newDevice).deviceID);
                obj.clBuilt=false;
                obj.clInitialized=false;
            end
        end
        
        %build the OpenCL program with selected precision, stepper, prob
        function buildCL(obj)
            obj.cppmethod('buildcl');
            obj.clBuilt=true;
        end
        
        %initialize builds the program and sets data needed to run
        %simulation in one call
        function initialize(obj, tspan, X0, P, sp)
            if ~exist('tspan','var') %no input args: use stored values
                tspan=obj.tspan; X0=obj.X0; P=obj.P; sp=obj.sp;
            end
            obj.cppmethod('initialize', tspan, X0(:), P(:), sp);
            obj.tspan=tspan;
            obj.X0=X0;
            obj.P=P;
            obj.sp=sp;
            obj.nPts=numel(X0)/obj.prob.nVar;
            obj.clInitialized=true;
        end
        
        function setNPts(obj, newNPts)
            if ~exist('newNPts','var') %no input args: use stored values
                newNPts=obj.nPts;
            end
            obj.cppmethod('setnpts', newNPts);
            obj.nPts=newNPts;
        end
        
        %Set X0 and P together if trying to change nPts
        function setProblemData(obj, X0, P)
            if ~exist('X0','var') %no input args: use stored values
                X0=obj.X0; P=obj.P;
            end
            obj.cppmethod('setproblemdata', X0(:), P(:));
            obj.X0=X0;
            obj.P=P;
            obj.nPts=numel(X0)/obj.prob.nVar;
        end
        
        
        function seedRNG(obj, mySeed)
            if ~exist('mySeed','var')
                obj.cppmethod('seedrng');
            else
                obj.cppmethod('seedrng', mySeed);
            end
        end
        
        
        function settspan(obj, tspan)
            if ~exist('tspan','var') %no input args: use stored values
                tspan=obj.tspan;
            end
            obj.cppmethod('settspan', tspan);
            obj.tspan=tspan;
        end
        
        %nPts cannot change here
        function setX0(obj, X0)
            if ~exist('X0','var') %no input args: use stored values
                X0=obj.X0;
            end
            testnPts=numel(X0)/obj.prob.nVar;
            if testnPts==obj.nPts %this check is redundant: c++ does it too
                obj.cppmethod('setx0', X0(:));
                obj.X0=X0;
            else
                error('Size of X0 is incorrect');
            end
        end
        
        %nPts cannot change here
        function setP(obj, P)
            if ~exist('P','var') %no input args: use stored values
                P=obj.P;
            end
            testnPts=numel(P)/obj.prob.nPar;
            if testnPts==obj.nPts
                obj.cppmethod('setpars', P(:));
                obj.P=P;
            else
                error('Size of P is incorrect');
            end
        end
        
        function setSolverPars(obj, sp)
            if ~exist('sp','var') %no input args: use stored values
                sp=obj.sp;
            end
            obj.cppmethod('setsolverpars', sp);
            obj.sp=sp;
        end
        
        % Integration
        function Xf=transient(obj, tspan)
            if ~obj.clBuilt
                error('OpenCL program not built. run buildCL')
            end
            if exist('tspan','var')
                obj.settspan(tspan);
            end
            obj.cppmethod('transient');
            if nargout==1 %overload to also transfer Xf from device to host
                Xf=obj.getXf();
            end
        end
        
        function shiftTspan(obj)
            obj.cppmethod('shifttspan');
            obj.getTspan();
        end
        
        function X0=shiftX0(obj)
            obj.cppmethod('shiftx0'); %device-device transfer Xf to X0
            if nargout==1 %overload to also transfer X0 from device to host
                X0=obj.getX0();
            end
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
    
    
    %ui public methods
    methods (Access=public)
        
        function hpar=uisetpars(obj,parent)
            if ~exist('parent','var')||isempty(parent)
                parent=uifigure();
            end
            partable=struct2table(obj.prob.par);
            hpar=uitable(parent,'Data',partable);
            hpar.ColumnName = partable.Properties.VariableNames;
%             hpar.ColumnName = {'name'; 'value'; 'lb'; 'ub'};
            hpar.ColumnWidth = {'auto', 'fit', 'fit', 'fit'};
            hpar.RowName = {};
            hpar.ColumnSortable = [false false false false];
            hpar.ColumnEditable = [false true true true];
%             hpar.Position = position;
            hpar.CellEditCallback=@updatePars;
            
            function updatePars(src,event)
                row=event.Indices(1);
                col=event.Indices(2);
                fieldname=hpar.DisplayData.Properties.VariableNames{col};
                if isnumeric(event.NewData) && event.NewData>0
                    obj.prob.par(row).(fieldname)=event.NewData;
                    if fieldname=="value"
                        obj.prob.p0(row)=event.NewData;
                    end
                end
            end
        end
        
        function hic=uisetIC(obj,parent)
            if ~exist('parent','var')||isempty(parent)
                parent=uifigure();
            end
            fullvartable=struct2table(obj.prob.var);
            vartable=fullvartable(:,1:4);
            hic=uitable(parent,'Data',vartable);
            hic.ColumnName = vartable.Properties.VariableNames;
%             hic.ColumnName = {'name'; 'value'; 'lb'; 'ub'};
            hic.ColumnWidth = {'auto', 'fit', 'fit', 'fit'};
            hic.RowName = {};
            hic.ColumnSortable = [false false false false];
            hic.ColumnEditable = [false true true true];
%             hpar.Position = position;
            hic.CellEditCallback=@updateVars;
            
            function updateVars(src,event)
                row=event.Indices(1);
                col=event.Indices(2);
                fieldname=hic.DisplayData.Properties.VariableNames{col};
                if isnumeric(event.NewData) && event.NewData>0
                    obj.prob.var(row).(fieldname)=event.NewData;
                    if fieldname=="value"
                        obj.prob.x0(row)=event.NewData;
                    end
                end
            end
        end
        
        function hsp=uisetSP(obj,parent)
            if ~exist('parent','var')||isempty(parent)
                parent=uifigure();
            end
            %TODO: convert XPP opts to name/value struct...
            %TODO: subset to relevant params (adaptive vs not)
            sptable=table;
            sptable.name=fieldnames(obj.sp);
            sptable.value=cell2mat(struct2cell(obj.sp));
            hsp=uitable(parent,'Data',sptable);
            hsp.ColumnName = sptable.Properties.VariableNames;
%             hpar.ColumnName = {'name'; 'value'};
            hsp.ColumnWidth = {'auto', 'fit'};
            hsp.RowName = {};
            hsp.ColumnSortable = [false false];
            hsp.ColumnEditable = [false true];
%             hsp.Position = position;
            hsp.CellEditCallback=@updateSP;
            
            function updateSP(src,event)
                thisfield=hsp.DisplayData.name{event.Indices(1)};
                if isnumeric(event.NewData) && event.NewData>0
                    obj.sp.(thisfield)=event.NewData;
                end
            end
        end
    end
    
    %ui private methods
    methods (Access=private)
    end
    
    %static helper methods
    methods (Static=true)
        
        function newprob = purifyProblemStruct(prob)
            newprob = struct;
            newprob.clRHSfilename = char(prob.clRHSfilename);
            newprob.nVar = prob.nVar;
            newprob.nPar = prob.nPar;
            newprob.nAux = prob.nAux;
            newprob.nWiener = prob.nWiener;
            newprob.varNames = prob.varNames;
            newprob.parNames = prob.parNames;
            newprob.auxNames = prob.auxNames;
        end

        function sp=defaultSolverParams()
            sp.dt=.1;
            sp.dtmax=100.00;
            sp.abstol=1e-6;
            sp.reltol=1e-4;
            sp.max_steps=1000000;
            sp.max_store=100000; %allocated number of timepoints: min( (tf-t0)/(dt*nout)+1 , sp.max_store)
            sp.nout=1;
%             sp.storevars=[]; %empty => all, otherwise specify list of var indices to store (for trajectories, or max/min/mean for features)
        end
        
        function selectedDevice=autoselectDevice(devices)
            selectedDevice=[];
%             if  isempty(selectedDevice)
%                 selectedDevice=find({devices(:).type}=="Accelerator",1,'first');
%             end
            if  isempty(selectedDevice)
                selectedDevice=find({devices(:).type}=="GPU",1,'first');
            end
            if  isempty(selectedDevice)
                selectedDevice=find({devices(:).type}=="CPU",1,'first');
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

