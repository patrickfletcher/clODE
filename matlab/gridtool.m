classdef gridtool < handle %matlab.mixin.SetGet
    %GRIDTOOL is a GUI for exploring the dependence of solutions of ODEs on
    %parameters, two-parameter planes at a time.  The ODEs are solved in
    %parallel using OpenCL via MEX interface to the c++ clODE library. The
    %interface is modeled after Bard Ermentrout's XPPAUT: ODE files that
    %run in XPP are parsed and converted to the required OpenCL code
    %automatically, and the interface shares some of the same keyboard
    %shortcuts.
    
    %properties representing clODE and the initial value problem
    properties (Access = public) %TODO: make private later
        
        odefile %ode system definition (an xpp ODE file)
        
        devices
        grid_device
        clo_g %clODE grid feature solver
        traj_device
        clo_t %clODE trajectory solver
        
        prob %read-only copy of the parseODE output
   
        trajvars %table of possible trajectory varibles (for easier lookup)
        %.name, type, ix
            
        nGrid=[32,32] %nPts in x, y dims
        nZ=10
        nPts=1024
        P0 %base parameter matrix with p0 replicated nPts times [nPts, nPar]
        X0 %base ic matrix with x0 replicated nPts times
        
        nClick=3
        
        %mechanism to cache F(:,:,z) values - for rapidly scanning z post compute
        
    end
    
    %data needed to display graphs: set-observable, so that whenever these
    %change, triggers update to relevant figure.
    properties %(setObservable = true)
        %should gridvars go here? update gridvars-->makeGridData & makeTrajData
        grid %struct specifying grid x, y, z info (name, N, del, vals)
        gridvars %working table of possible grid vars (par+var) with values
        gridSol=struct('x',[],'y',[],'F',[]);
        trajP0Sol=struct('t',[],'x',[],'aux',[]); %struct for p0 trajectory data
        trajClickSol=struct('t',[],'x',[],'aux',[]); %array of nClick structs for pQuery trajectory data
    end
    
    %listeners
    properties (Access = private)
        listenerGridData
        listenerTrajP0
        listenerTrajClick
    end
    
    % Properties that correspond to UI components
    properties (Access = public)
        
        %UI state
        gridclBuilt=false %flag to indicate if OpenCL program is built
        trajclBuilt=false
        gridclInitialized=false %flag to indicate if ready for integration
        trajclInitialized=false
        
        % change in grid, par.value, x0.value, sp, op => solution not current
        gridIsCurrent=false  %flag to indicate need for new solution
        trajP0SolIsCurrent=false
        trajClickIsCurrent=false
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %grid figure
        figGrid                 
        axGrid                  matlab.graphics.axis.Axes
        imGrid                  matlab.graphics.primitive.Image
        gridCBar
        markerP0 %marker (line object) specifying (x,y) for p0
        markerPquery %markers (line objects) specifying (x,y) for n query trajectories

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %figure with p0 trajectory
        figTrajP0               
        axTrajP0                matlab.graphics.axis.Axes
        lineTrajP0              matlab.graphics.chart.primitive.Line  %could be array of two if yyaxis
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %click trajectory figure with n trajectories
        figTrajClick  
        tilesTrajClick
        axTrajClick             matlab.graphics.axis.Axes
        lineTrajClick           matlab.graphics.chart.primitive.Line  %could be array of 2*nClick if yyaxis
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %control figure
        figControl              matlab.ui.Figure
        TabGroup                matlab.ui.container.TabGroup
        
        NumericsTab             matlab.ui.container.Tab
        
        SelectODEfileButton     matlab.ui.control.Button
        currentFileLabel        matlab.ui.control.Label
        
        gridDeviceDropDownLabel   matlab.ui.control.Label
        gridDeviceDropDown        matlab.ui.control.DropDown
        gridPrecisionCheckBox     matlab.ui.control.CheckBox
        gridStepperDropDownLabel  matlab.ui.control.Label
        gridStepperDropDown       matlab.ui.control.DropDown
        gridObserverDropDownLabel matlab.ui.control.Label
        gridObserverDropDown      matlab.ui.control.DropDown
        gridCLIsBuiltButton       matlab.ui.control.StateButton
        
        trajDeviceDropDownLabel   matlab.ui.control.Label
        trajStepperDropDownLabel  matlab.ui.control.Label
        trajStepperDropDown       matlab.ui.control.DropDown
        trajPrecisionCheckBox     matlab.ui.control.CheckBox
        trajDeviceDropDown        matlab.ui.control.DropDown
        trajCLIsBuiltButton       matlab.ui.control.StateButton
        
        gridSolverParLabel      matlab.ui.control.Label
        gridSolverParTable      matlab.ui.control.Table
        trajSolverParTable      matlab.ui.control.Table
        trajSolverParLabel      matlab.ui.control.Label
        gridObserverParTable    matlab.ui.control.Table
        gridObserverParLabel    matlab.ui.control.Label
        
        GridTab                 matlab.ui.container.Tab
        featureDropDownLabel    matlab.ui.control.Label
        featureDropDown         matlab.ui.control.DropDown
        fscaleEditFieldLabel    matlab.ui.control.Label
        fscaleEditField         matlab.ui.control.NumericEditField
        
        gridTable               matlab.ui.control.Table
        tspanTable              matlab.ui.control.Table
        parLabel                matlab.ui.control.Label
        parTable                matlab.ui.control.Table
        icLabel                 matlab.ui.control.Label
        icTable                 matlab.ui.control.Table
        
        parDefaultButton        matlab.ui.control.Button
        icDefaultButton         matlab.ui.control.Button
        icRandomButton          matlab.ui.control.Button
        trajectoriesLabel       matlab.ui.control.Label
        clicksEditFieldLabel    matlab.ui.control.Label
        clicksEditField         matlab.ui.control.NumericEditField
        xDropDownLabel          matlab.ui.control.Label
        xDropDown               matlab.ui.control.DropDown
        yDropDownLabel          matlab.ui.control.Label
        yDropDown               matlab.ui.control.DropDown
        zDropDownLabel          matlab.ui.control.Label
        zDropDown               matlab.ui.control.DropDown
        linkAxesButton          matlab.ui.control.StateButton
        threeDButton            matlab.ui.control.StateButton
        
%         HelpTab                 matlab.ui.container.Tab
    end

    % Callbacks that handle component events
    methods (Access = private)
        
        % Close request function: figControl
        function figControlCloseRequest(app, src, event)
            close(app.figGrid)
            close(app.figTrajP0)
            close(app.figTrajClick)
            closereq
%             delete(app);
%             clear gridtool
        end
        
        % Button pushed function: SelectODEfileButton
        function SelectODEfileButtonPushed(app, src, event)
            % choose ODE: trigger ode2cl, populate all components with data 
            [name,path]=uigetfile('.ode','Select an ODE file');
            if ~ischar(name)
                disp('File selection canceled, quitting...')
                return
            end
            app.odefile=fullfile(path,name);
            app.processNewODEfile(app.odefile);
        end

        % Value changed function: gridDeviceDropDown
        function gridDeviceDropDownValueChanged(app, src, event)
            value = app.gridDeviceDropDown.Value;
            app.clo_g.setOpenCL(value);
            app.gridclInitialized=false;
        end

        % Value changed function: trajDeviceDropDown
        function trajDeviceDropDownValueChanged(app, src, event)
            value = app.trajDeviceDropDown.Value;
            app.clo_t.setOpenCL(value); 
            app.trajclInitialized=false;           
        end

        % Value changed function: gridPrecisionCheckBox
        function gridPrecisionCheckBoxValueChanged(app, src, event)
            useDouble = app.gridPrecisionCheckBox.Value;
            precision='single';
            if useDouble, precision='double'; end
            app.clo_g.setPrecision(precision);
            app.gridclInitialized=false;
        end

        % Value changed function: trajPrecisionCheckBox
        function trajPrecisionCheckBoxValueChanged(app, src, event)
            useDouble = app.trajPrecisionCheckBox.Value;
            precision='single';
            if useDouble, precision='double'; end
            app.clo_t.setPrecision(precision);
            app.trajclInitialized=false;
        end

        % Value changed function: gridStepperDropDown
        function gridStepperDropDownValueChanged(app, src, event)
            stepper = app.gridStepperDropDown.Value;
            app.clo_g.setStepper(stepper);
            app.gridclInitialized=false;
        end

        % Value changed function: trajStepperDropDown
        function trajStepperDropDownValueChanged(app, src, event)
            stepper = app.trajStepperDropDown.Value;
            app.clo_t.setStepper(stepper);
            app.trajclInitialized=false;
        end
        
        % Value changed function: gridObserverDropDown
        function gridObserverDropDownValueChanged(app, src, event)
            observer = app.gridObserverDropDown.Value;
            app.clo_g.setObserver(observer);
            app.gridclInitialized=false;
        end

        % Value changed function: gridCLIsBuiltButton (doesn't execute for
        % programmatic changes)
        function gridCLIsBuiltButtonValueChanged(app, src, event)
            isBuilt = app.gridCLIsBuiltButton.Value;
            if ~isBuilt
                app.buildCL('grid');
            end
        end

    	% Value changed function: trajCLIsBuiltButton (doesn't execute for
        % programmatic changes)
        function trajCLIsBuiltButtonValueChanged(app, src, event)
            isBuilt = app.trajCLIsBuiltButton.Value;
            if ~isBuilt
                app.buildCL('traj');
            end
        end

        % helper functions for build/build needed
        function buildCL(app,which_clo)
            switch which_clo
                case 'grid'
                    app.clo_g.buildCL();
                    app.gridCLIsBuiltButton.Value = 1;
                    app.gridCLIsBuiltButton.Text = 'CL ready';
                    app.gridCLIsBuiltButton.BackgroundColor=[0.7,1,0.7];
                case 'traj'
                    app.clo_t.buildCL();
                    app.trajCLIsBuiltButton.Value = 1;
                    app.trajCLIsBuiltButton.Text = 'CL ready';
                    app.trajCLIsBuiltButton.BackgroundColor=[0.7,1,0.7];
            end
        end
        
        function setBuildNeeded(app,which_clo)
            switch which_clo
                case 'grid'
                    app.gridCLIsBuiltButton.Value = 0;
                    app.gridCLIsBuiltButton.Text = 'Build CL';
                    app.gridCLIsBuiltButton.BackgroundColor=[1,0.7,0.7];
                case 'traj'
                    app.trajCLIsBuiltButton.Value = 0;
                    app.trajCLIsBuiltButton.Text = 'Build CL';
                    app.trajCLIsBuiltButton.BackgroundColor=[1,0.7,0.7];
            end
        end
        
        
        % Cell edit callback: gridSolverParTable
        function gridSolverParTableCellEdit(app, src, event)
            indices = event.Indices;
            newData = event.NewData;
            thisfield=src.DisplayData.Properties.RowNames{indices(1)};
            if isnumeric(newData) && newData>0
                sp=app.clo_g.sp;
                sp.(thisfield)=newData;
                app.clo_g.setSolverPars(sp); %do the device transfer (could clean this up if clODE had set.sp)
            end
        end

        % Cell edit callback: trajSolverParTable
        function trajSolverParTableCellEdit(app, src, event)
            indices = event.Indices;
            newData = event.NewData;
            thisfield=src.DisplayData.Properties.RowNames{indices(1)};
            if isnumeric(newData) && newData>0
                sp=app.clo_t.sp;
                sp.(thisfield)=newData;
                app.clo_t.setSolverPars(sp);
            end
        end

        % Cell edit callback: gridObserverParTable
        function gridObserverParTableCellEdit(app, src, event)
            indices = event.Indices;
            newData = event.NewData;
            thisfield=src.DisplayData.Properties.RowNames{indices(1)};
            if isnumeric(newData) && newData>0
                app.clo_g.op.(thisfield)=newData;
                op=app.clo_g.op;
                op.(thisfield)=newData;
                app.clo_g.setObserverPars(op);
            end
        end

        % Value changed function: featureDropDown
        function featureDropDownValueChanged(app, src, event)
            newFeature = app.featureDropDown.Value;
            %trigger update grid
        end

        % Value changed function: fscaleEditField
        function fscaleEditFieldValueChanged(app, src, event)
            newFscale = app.fscaleEditField.Value;
            %trigger update grid
        end

        % Value changed function: t0EditField
        function t0EditFieldValueChanged(app, src, event)
            t0 = app.t0EditField.Value;
            if t0<app.clo_g.tspan(2)
                tspan(1)=t0;
                app.clo_g.settspan(tspan);
                app.clo_t.settspan(tspan);
            end
        end

        % Value changed function: tfEditField
        function tfEditFieldValueChanged(app, src, event)
            tf = app.tfEditField.Value;
            if tf>app.clo_g.tspan(1)
                tspan(2)=tf;
                app.clo_g.settspan(tspan);
                app.clo_t.settspan(tspan);
            end
        end

        % Cell edit callback: gridTable
        function gridTableCellEdit(app, src, event)
            %update clo_g.P, clo_g.X0. Link any lb/ub changes to par/var. support
            %changing del or N.
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Cell edit callback: parTable
        function parTableCellEdit(app, src, event)
            %update clo_g.P. Link any lb/ub changes to gridTable.
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Cell edit callback: icTable
        function icTableCellEdit(app, src, event)
            %update clo_g.X0. Link any lb/ub changes to gridTable.
            indices = event.Indices;
            newData = event.NewData;
        end

        % Button pushed function: parDefaultButton
        function parDefaultButtonPushed(app, src, event)
            %get the original default values for par.val/lb/ub from
            %app.prob
        end

        % Button pushed function: icDefaultButton
        function icDefaultButtonPushed(app, src, event)
            %get the original default values for var.val/lb/ub from
            %app.prob
        end

        % Button pushed function: icRandomButton
        function icRandomButtonPushed(app, src, event)
            %randomize between lb/ub for ic: update X0
        end

        % Value changed function: xDropDown
        function xDropDownValueChanged(app, src, event)
            newX = app.xDropDown.Value;
            
        end

        % Value changed function: yDropDown
        function yDropDownValueChanged(app, src, event)
            newY = app.yDropDown.Value;
            
        end

        % Value changed function: zDropDown
        function zDropDownValueChanged(app, src, event)
            newZ = app.zDropDown.Value;
            
        end

        % Value changed function: clicksEditField
        function clicksEditFieldValueChanged(app, src, event)
            newNClick = app.clicksEditField.Value;
            app.nClick=newNClick;
            
        end

        % Value changed function: linkAxesButton
        function linkAxesButtonValueChanged(app, src, event)
            doLinkAxes = app.linkAxesButton.Value;
        end

        % Value changed function: threeDButton
        function threeDButtonValueChanged(app, src, event)
            do3D = app.threeDButton.Value;
%             app.threeDButton
        end
    end
    
    % internal functions (clODE control & ui interaction)
    methods (Access = private)
        
        function gridKeyPress(app)
        end
        
        function trajKeyPress(app)
        end
        
        function makeGridData(app)
            newP=app.P0;
            newX0=app.X0;
            
            vals={};
            for i=1:3
                vals{i,1}=linspace(app.grid.lb(1),app.grid.ub(1),app.nGrid(1));
            end
            app.grid.vals=vals;
            [X,Y]=meshgrid(x,y);
            
            if app.grid{'x','type'}=="par"
                newP(:,app.grid{'x','ix'})=X(:);
            else
                newX0(:,app.grid{'x','ix'})=X(:);
            end
            if app.grid{'y','type'}=="par"
                newP(:,app.grid{'y','ix'})=Y(:);
            else
                newX0(:,app.grid{'y','ix'})=Y(:);
            end
        end
        
        function makeTrajData(app, clickCoords)
        end
        
        function integrateGrid(app)
%             app.clo_g.initialize(); %send all data to GPU.   NO: assume it is stored prior to this

        end
        
        function integrateTraj(app)
        end
        
        function updateGridPlot(app)
        end
        
        function updateTrajP0Plot(app)
        end
        
        function updateTrajClickPlot(app)
        end
        
        function processNewODEfile(app, odefile, firstTime)
            
            if ~exist('firstTime','var')
                firstTime=false;
            end
            
            app.odefile=odefile;
            [~,app.prob]=ode2cl(odefile);
            
            if app.prob.nPar+app.prob.nVar<2
                error('Not enough parameters and variables for a 2D grid!')
            end
            
            app.currentFileLabel.Text=app.prob.name+".ode";
            
            %If first time, set up clODE objects using their default parameters
            if firstTime
                app.clo_g=clODEfeatures(app.prob, [], app.grid_device);
                app.clo_t=clODEtrajectory(app.prob, [], app.traj_device);
                
                app.gridObserverDropDown.Items=app.clo_g.observerNames;
                app.gridObserverDropDown.Value=app.clo_g.observer;
                app.gridStepperDropDown.Items=app.clo_g.getAvailableSteppers();
                app.gridStepperDropDown.Value=app.clo_g.stepper;
                app.trajStepperDropDown.Items=app.clo_t.getAvailableSteppers();
                app.trajStepperDropDown.Value=app.clo_t.stepper;
                
                defaultsp=clODE.defaultSolverParams(); %app.clo_g.sp
                sptable=table('RowNames',fieldnames(defaultsp));
                sptable.value=cell2mat(struct2cell(defaultsp));
                app.gridSolverParTable.Data=sptable;
                app.trajSolverParTable.Data=sptable;
%                 app.clo_g.setSolverPars(sptable); %this is done later in init
%                 app.clo_t.setSolverPars(sptable); 
            
                defaultop=clODEfeatures.defaultObserverParams(); %app.clo_g.op
                optable=table('RowNames',fieldnames(defaultop));
                optable.value=cell2mat(struct2cell(defaultop));
                app.gridObserverParTable.Data=optable;
                app.gridObserverParTable.RowName=fieldnames(defaultop);
%                 app.clo_g.setObserverPars(defaultop); %this is done later in init
                
                app.featureDropDown.Items=app.clo_g.featureNames();
                
            else %to avoid reverting to default sp and op:
                app.clo_g.setNewProblem(app.prob);
                app.clo_t.setNewProblem(app.prob);
            end
            app.buildCL('grid');
            app.buildCL('traj');
            
            app.gridclInitialized=false;
            app.trajclInitialized=false;
            
            %solver opts from ODE file
            dt=app.prob.opt.dt;
            %dtmax, abstol/reltol, nout, maxstore?
            app.gridSolverParTable.Data{'dt','value'}=dt;
            app.trajSolverParTable.Data{'dt','value'}=dt;
            
            %tspan
            t0=app.prob.opt.t0;
            tf=app.prob.opt.total;
            app.clo_g.settspan([t0,tf]); %device transfer - defer to init
            app.clo_t.settspan([t0,tf]); %device transfer
            
            tspan=table;
            tspan.t0=[t0;t0];
            tspan.tf=[tf;tf];
            tspan.trans=[tf;tf];
            app.tspanTable.Data=tspan;
            
            %valid grid variables
            icNames=strcat(app.prob.varNames(:),'0');
            app.gridvars=table;
            app.gridvars.name=[app.prob.parNames(:);icNames(:)];
            app.gridvars.value=[app.prob.p0(:);app.prob.x0(:)];
            app.gridvars.lb=[[app.prob.par.lb]';[app.prob.var.lb]'];
            app.gridvars.ub=[[app.prob.par.ub]';[app.prob.var.ub]'];
            app.gridvars.type=categorical([ones(app.prob.nPar,1);2*ones(app.prob.nVar,1)],[1,2],{'par','ic'});
            app.gridvars.ix=[(1:app.prob.nPar)';(1:app.prob.nVar)'];
            app.gridvars.Properties.RowNames=app.gridvars.name;
            
            % trajectory variables
            app.trajvars=table;
            app.trajvars.name=[app.prob.varNames(:);app.prob.auxNames(:)];
            app.trajvars.type=categorical([ones(app.prob.nVar,1);2*ones(app.prob.nAux,1)],[1,2],{'var','aux'});
            app.trajvars.ix=[(1:app.prob.nVar)';(1:app.prob.nAux)'];
            app.trajvars.Properties.RowNames=app.trajvars.name;
            
            %set up initial default gridTable
            app.grid=app.gridvars.name(1:3,:);
            app.grid.Properties.RowNames={'x','y','z'};
            app.grid.N=[app.nGrid(:);app.nZ];
            for i=1:3 
                var=app.grid.name(i);
                lb{i,1}=app.gridvars{var,'lb'};
                ub{i,1}=app.gridvars{var,'ub'};
                app.grid.del(i,1)=(app.grid.ub(i)-app.grid.lb(i))/(app.grid.N(i)-1); %compute dels
            end
            app.gridTable.Data=table2cell(app.grid(:,{'name','N','del'}));
            app.gridTable.ColumnFormat={app.gridvars.name',[],[]};
            xlabel(app.axGrid,app.grid.name(1));
            ylabel(app.axGrid,app.grid.name(2));
            
            %parameters
            app.parTable.Data=app.gridvars(app.gridvars.type=="par",{'name','value','lb','ub'});
            
            %ic table
            app.icTable.Data=app.gridvars(app.gridvars.type=="ic",{'name','value','lb','ub'});
            
            %trajectory plot variables
            app.xDropDown.Items=[{'t'};app.trajvars.name];
            app.yDropDown.Items=app.trajvars.name;
            app.zDropDown.Items=[{'none'};app.trajvars.name];
            
            xlabel(app.axTrajP0,app.xDropDown.Value);
            ylabel(app.axTrajP0,app.yDropDown.Value);
            xlabel(app.axTrajClick(end),app.xDropDown.Value);
            ylabel(app.axTrajClick(end),app.yDropDown.Value);
            
            %make the actual grid data for simulation
            app.P0=repmat(app.gridvars.value(app.gridvars.type=="par")',app.nPts,1);
            app.X0=repmat(app.gridvars.value(app.gridvars.type=="var")',app.nPts,1);
            app.makeGridData();
            app.makeTrajData(); %p0: no clickCoords
            
            figure(app.figControl)
        end
        
    end

    
    methods (Access = public)
        %constructor
        function app=gridtool(odefile)%, precision, selectedDevice, stepper, observer
            
            app.createUIComponents();
            
            %OpenCL device list
            app.devices=queryOpenCL();
            
            %auto select first GPU for grid, fastest clock for traj
            app.grid_device=find({app.devices(:).type}=="GPU",1,'first');
            [~,app.traj_device]=max([app.devices(:).maxClock]);
            
            app.gridDeviceDropDown.Items={app.devices(:).name};
            app.gridDeviceDropDown.Value=app.devices(app.grid_device).name;
            app.trajDeviceDropDown.Items={app.devices(:).name};
            app.trajDeviceDropDown.Value=app.devices(app.traj_device).name;
            
            if exist('odefile','var')&&~isempty(odefile)
                app.processNewODEfile(odefile, 1)
            end
            
%             if nargout==0
%                 clear app
%             end

        end
        
    end
    
    
    % Component initialization
    methods (Access = private)

        % Create UI figures and components
        function createUIComponents(app)
            app.createControlFig();
            app.createGridFig();
            app.createTrajP0Fig();
            app.createTrajClickFig();
        end
        
        function createControlFig(app)
            %TODO: normalized position units, relative to screen size?
        
            % Create figControl and hide until all components are created
            app.figControl = uifigure('Visible', 'off');
            app.figControl.Position = [50 450 300 600];
            app.figControl.Name = 'gridtool';
            app.figControl.CloseRequestFcn = @app.figControlCloseRequest;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.figControl);
            app.TabGroup.Position = [5 5 290 590];

            % Create NumericsTab
            app.NumericsTab = uitab(app.TabGroup);
            app.NumericsTab.Title = 'Numerics';

            % Create SelectODEfileButton
            app.SelectODEfileButton = uibutton(app.NumericsTab, 'push');
            app.SelectODEfileButton.ButtonPushedFcn = @app.SelectODEfileButtonPushed;
            app.SelectODEfileButton.Position = [13 530 115 22];
            app.SelectODEfileButton.Text = 'Select ODE file';

            % Create currentFileLabel
            app.currentFileLabel = uilabel(app.NumericsTab);
            app.currentFileLabel.Position = [137 530 138 22];
            app.currentFileLabel.Text = 'currentFile';

            % Create gridDeviceDropDownLabel
            app.gridDeviceDropDownLabel = uilabel(app.NumericsTab);
            app.gridDeviceDropDownLabel.Position = [246 496 40 22];
            app.gridDeviceDropDownLabel.Text = 'device';

            % Create gridDeviceDropDown
            app.gridDeviceDropDown = uidropdown(app.NumericsTab);
            app.gridDeviceDropDown.ValueChangedFcn = @app.gridDeviceDropDownValueChanged;
            app.gridDeviceDropDown.Position = [172 478 114 22];

            % Create gridPrecisionCheckBox
            app.gridPrecisionCheckBox = uicheckbox(app.NumericsTab);
            app.gridPrecisionCheckBox.ValueChangedFcn = @app.gridPrecisionCheckBoxValueChanged;
            app.gridPrecisionCheckBox.Text = 'double precision';
            app.gridPrecisionCheckBox.Position = [172 452 109 23];

            % Create gridStepperDropDownLabel
            app.gridStepperDropDownLabel = uilabel(app.NumericsTab);
            app.gridStepperDropDownLabel.HorizontalAlignment = 'center';
            app.gridStepperDropDownLabel.Position = [172 427 46 22];
            app.gridStepperDropDownLabel.Text = 'stepper';

            % Create gridStepperDropDown
            app.gridStepperDropDown = uidropdown(app.NumericsTab);
            app.gridStepperDropDown.ValueChangedFcn = @app.gridStepperDropDownValueChanged;
            app.gridStepperDropDown.Position = [217 427 69 22];

            % Create gridSolverParLabel
            app.gridSolverParLabel = uilabel(app.NumericsTab);
            app.gridSolverParLabel.Position = [13 498 128 22];
            app.gridSolverParLabel.Text = 'Grid solver parameters';

            % Create gridSolverParTable
            app.gridSolverParTable = uitable(app.NumericsTab);
            app.gridSolverParTable.ColumnName = {'value'};
            app.gridSolverParTable.ColumnWidth = {'fit'};
            app.gridSolverParTable.ColumnEditable = [true];
            app.gridSolverParTable.CellEditCallback = @app.gridSolverParTableCellEdit;
            app.gridSolverParTable.Position = [13 355 145 145];
            


            % Create gridObserverDropDownLabel
            app.gridObserverDropDownLabel = uilabel(app.NumericsTab);
            app.gridObserverDropDownLabel.HorizontalAlignment = 'center';
            app.gridObserverDropDownLabel.Position = [172 330 52 22];
            app.gridObserverDropDownLabel.Text = 'observer';

            % Create gridObserverDropDown
            app.gridObserverDropDown = uidropdown(app.NumericsTab);
            app.gridObserverDropDown.ValueChangedFcn = @app.gridObserverDropDownValueChanged;
            app.gridObserverDropDown.Position = [175 309 107 22];
            
            % Create gridObserverParLabel
            app.gridObserverParLabel = uilabel(app.NumericsTab);
            app.gridObserverParLabel.Position = [13 330 119 22];
            app.gridObserverParLabel.Text = 'Observer parameters';

            % Create gridObserverParTable
            app.gridObserverParTable = uitable(app.NumericsTab);
            app.gridObserverParTable.ColumnName = {'value'};
            app.gridObserverParTable.ColumnWidth = {'fit'};
            app.gridObserverParTable.ColumnEditable = [true];
            app.gridObserverParTable.CellEditCallback = @app.gridObserverParTableCellEdit;
            app.gridObserverParTable.Position = [13 186 150 145];
            
            % Create gridCLIsBuiltButton
            app.gridCLIsBuiltButton = uibutton(app.NumericsTab, 'state');
            app.gridCLIsBuiltButton.ValueChangedFcn = @app.gridCLIsBuiltButtonValueChanged;
            app.gridCLIsBuiltButton.Text = 'Build CL';
            app.gridCLIsBuiltButton.Value=0;
            app.gridCLIsBuiltButton.BackgroundColor=[1,.7,.7];
            app.gridCLIsBuiltButton.Position = [217 186 69 22];
            

            % Create trajDeviceDropDownLabel
            app.trajDeviceDropDownLabel = uilabel(app.NumericsTab);
            app.trajDeviceDropDownLabel.Position = [246 161 40 22];
            app.trajDeviceDropDownLabel.Text = 'device';

            % Create trajDeviceDropDown
            app.trajDeviceDropDown = uidropdown(app.NumericsTab);
            app.trajDeviceDropDown.ValueChangedFcn = @app.trajDeviceDropDownValueChanged;
            app.trajDeviceDropDown.Position = [172 140 114 22];

            % Create trajPrecisionCheckBox
            app.trajPrecisionCheckBox = uicheckbox(app.NumericsTab);
            app.trajPrecisionCheckBox.ValueChangedFcn = @app.trajPrecisionCheckBoxValueChanged;
            app.trajPrecisionCheckBox.Text = 'double precision';
            app.trajPrecisionCheckBox.Position = [172 112 109 23];

            % Create trajStepperDropDownLabel
            app.trajStepperDropDownLabel = uilabel(app.NumericsTab);
            app.trajStepperDropDownLabel.HorizontalAlignment = 'center';
            app.trajStepperDropDownLabel.Position = [172 86 46 22];
            app.trajStepperDropDownLabel.Text = 'stepper';

            % Create trajStepperDropDown
            app.trajStepperDropDown = uidropdown(app.NumericsTab);
            app.trajStepperDropDown.ValueChangedFcn = @app.trajStepperDropDownValueChanged;
            app.trajStepperDropDown.Position = [217 86 69 22];

            % Create trajSolverParLabel
            app.trajSolverParLabel = uilabel(app.NumericsTab);
            app.trajSolverParLabel.Position = [13 161 158 22];
            app.trajSolverParLabel.Text = 'Trajectory solver parameters';
            
            % Create trajSolverParTable
            app.trajSolverParTable = uitable(app.NumericsTab);
            app.trajSolverParTable.ColumnName = {'value'};
            app.trajSolverParTable.ColumnWidth = {'fit'};
            app.trajSolverParTable.ColumnEditable = [true];
            app.trajSolverParTable.CellEditCallback = @app.trajSolverParTableCellEdit;
            app.trajSolverParTable.Position = [15 15 145 145];
            
            % Create trajCLIsBuiltButton
            app.trajCLIsBuiltButton = uibutton(app.NumericsTab, 'state');
            app.trajCLIsBuiltButton.ValueChangedFcn = @app.trajCLIsBuiltButtonValueChanged;
            app.trajCLIsBuiltButton.Text = 'Build CL';
            app.trajCLIsBuiltButton.Value=0;
            app.trajCLIsBuiltButton.BackgroundColor=[1,.7,.7];
            app.trajCLIsBuiltButton.Position = [217 15 69 22];

            % Create GridTab
            app.GridTab = uitab(app.TabGroup);
            app.GridTab.Title = 'Grid';
            
            % Create featureDropDownLabel
            app.featureDropDownLabel = uilabel(app.GridTab);
            app.featureDropDownLabel.Position = [15 530 43 22];
            app.featureDropDownLabel.Text = 'feature';

            % Create featureDropDown
            app.featureDropDown = uidropdown(app.GridTab);
            app.featureDropDown.ValueChangedFcn = @app.featureDropDownValueChanged;
            app.featureDropDown.Position = [61 530 106 22];

            % Create fscaleEditFieldLabel
            app.fscaleEditFieldLabel = uilabel(app.GridTab);
            app.fscaleEditFieldLabel.HorizontalAlignment = 'center';
            app.fscaleEditFieldLabel.Position = [176 530 43 22];
            app.fscaleEditFieldLabel.Text = 'fscale';

            % Create fscaleEditField
            app.fscaleEditField = uieditfield(app.GridTab, 'numeric');
            app.fscaleEditField.ValueChangedFcn = @app.fscaleEditFieldValueChanged;
            app.fscaleEditField.Position = [220 532 55 18];
            app.fscaleEditField.Value=1;
            

            % Create gridTable
            app.gridTable = uitable(app.GridTab);
            app.gridTable.RowName={'x';'y';'z'};
            app.gridTable.ColumnName={'name','N','del'};
            app.gridTable.ColumnWidth = {'fit', 'fit', 'fit'};
            app.gridTable.ColumnEditable = [true true false];
            app.gridTable.CellEditCallback = @app.gridTableCellEdit;
            app.gridTable.Position = [15 430 150 100];
            
            % Create tspanTable
            app.tspanTable = uitable(app.GridTab);
            app.tspanTable.ColumnName = {'t0','tf','trans'};
            app.tspanTable.RowName = {'grid','trajectory'};
            app.tspanTable.ColumnWidth = {'fit', 'fit', 'fit'};
            app.tspanTable.ColumnEditable = [true true true];
            app.tspanTable.CellEditCallback = @app.gridTableCellEdit;
            app.tspanTable.Position = [170 430 125 100];
            

            % Create parLabel
            app.parLabel = uilabel(app.GridTab);
            app.parLabel.Position = [15 403 69 23];
            app.parLabel.Text = 'Parameters';

            % Create parTable
            app.parTable = uitable(app.GridTab);
            app.parTable.RowName={};
            app.parTable.ColumnName = {'name', 'value', 'lb', 'ub'};
            app.parTable.ColumnWidth = {'fit', 'fit', 'fit', 'fit'};
            app.parTable.ColumnSortable = [true false false false];
            app.parTable.ColumnEditable = [true true true true];
            app.parTable.CellEditCallback = @app.parTableCellEdit;
            app.parTable.Position = [15 245 205 160];

            % Create parDefaultButton
            app.parDefaultButton = uibutton(app.GridTab, 'push');
            app.parDefaultButton.ButtonPushedFcn = @app.parDefaultButtonPushed;
            app.parDefaultButton.Position = [229 382 52 22];
            app.parDefaultButton.Text = 'default';
            

            % Create icLabel
            app.icLabel = uilabel(app.GridTab);
            app.icLabel.Position = [15 219 91 23];
            app.icLabel.Text = 'Initial conditions';

            % Create icTable
            app.icTable = uitable(app.GridTab);
            app.icTable.RowName={};
            app.icTable.ColumnName = {'name', 'value', 'lb', 'ub'};
            app.icTable.ColumnWidth = {'fit', 'fit', 'fit', 'fit'};
            app.icTable.ColumnSortable = [true false false false];
            app.icTable.ColumnEditable = [true true true true];
            app.icTable.CellEditCallback = @app.icTableCellEdit;
            app.icTable.Position = [15 70 205 150];

            
            % Create icDefaultButton
            app.icDefaultButton = uibutton(app.GridTab, 'push');
            app.icDefaultButton.ButtonPushedFcn = @app.icDefaultButtonPushed;
            app.icDefaultButton.Position = [229 198 52 22];
            app.icDefaultButton.Text = 'default';

            % Create icRandomButton
            app.icRandomButton = uibutton(app.GridTab, 'push');
            app.icRandomButton.ButtonPushedFcn = @app.icRandomButtonPushed;
            app.icRandomButton.Position = [229 172 52 22];
            app.icRandomButton.Text = 'random';
            

            % Create trajectoriesLabel
            app.trajectoriesLabel = uilabel(app.GridTab);
            app.trajectoriesLabel.Position = [15 34 68 23];
            app.trajectoriesLabel.Text = 'Trajectories';

            % Create xDropDownLabel
            app.xDropDownLabel = uilabel(app.GridTab);
            app.xDropDownLabel.HorizontalAlignment = 'right';
            app.xDropDownLabel.Position = [15 7 10 22];
            app.xDropDownLabel.Text = 'x';
            
            % Create xDropDown
            app.xDropDown = uidropdown(app.GridTab);
            app.xDropDown.ValueChangedFcn = @app.xDropDownValueChanged;
            app.xDropDown.Position = [32 7 52 21];

            % Create yDropDownLabel
            app.yDropDownLabel = uilabel(app.GridTab);
            app.yDropDownLabel.HorizontalAlignment = 'right';
            app.yDropDownLabel.Position = [83 7 10 22];
            app.yDropDownLabel.Text = 'y';

            % Create yDropDown
            app.yDropDown = uidropdown(app.GridTab);
            app.yDropDown.ValueChangedFcn = @app.yDropDownValueChanged;
            app.yDropDown.Position = [100 7 52 21];

            % Create zDropDownLabel
            app.zDropDownLabel = uilabel(app.GridTab);
            app.zDropDownLabel.HorizontalAlignment = 'right';
            app.zDropDownLabel.Position = [152 7 10 22];
            app.zDropDownLabel.Text = 'z';

            % Create zDropDown
            app.zDropDown = uidropdown(app.GridTab);
            app.zDropDown.ValueChangedFcn = @app.zDropDownValueChanged;
            app.zDropDown.Position = [169 7 52 21];

            % Create clicksEditFieldLabel
            app.clicksEditFieldLabel = uilabel(app.GridTab);
            app.clicksEditFieldLabel.HorizontalAlignment = 'right';
            app.clicksEditFieldLabel.Position = [146 35 44 21];
            app.clicksEditFieldLabel.Text = '# clicks';

            % Create clicksEditField
            app.clicksEditField = uieditfield(app.GridTab, 'numeric');
            app.clicksEditField.ValueChangedFcn = @app.clicksEditFieldValueChanged;
            app.clicksEditField.Position = [193 34 27 21];
            app.clicksEditField.Value=app.nClick;

            % Create linkAxesButton
            app.linkAxesButton = uibutton(app.GridTab, 'state');
            app.linkAxesButton.ValueChangedFcn = @app.linkAxesButtonValueChanged;
            app.linkAxesButton.Text = 'link axes';
            app.linkAxesButton.Position = [227 34 55 22];

            % Create threeDButton
            app.threeDButton = uibutton(app.GridTab, 'state');
            app.threeDButton.ValueChangedFcn = @app.threeDButtonValueChanged;
            app.threeDButton.Text = '3D';
            app.threeDButton.Position = [227 7 55 22];
            
            
%             app.HelpTab = uitab(app.TabGroup);
%             app.HelpTab.Title = 'Usage';
            
            
            app.figControl.Visible = 'on';
            
        end
            
        function createGridFig(app)
            % Create figGrid and hide until all components are created
            app.figGrid = figure('Visible', 'off');
            app.figGrid.Position = [375 450 600 600];
            app.figGrid.Name = 'grid';
            
            % grid axis
            t=tiledlayout(app.figGrid,1,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axGrid = nexttile(t);
            
            app.imGrid=imagesc(app.axGrid,0, 0, 0);
            app.imGrid.Visible='off';
            app.imGrid.HitTest='off';
            app.axGrid.YDir='normal';
            app.gridCBar=colorbar('northoutside');
%             app.gridCBar.Position

            xlabel(app.axGrid,'x')
            ylabel(app.axGrid,'y')
            
            app.figGrid.Visible = 'on';
        end
        
        function createTrajP0Fig(app)
            % Create figTraj and hide until all components are created
            app.figTrajP0 = figure('Visible', 'off');
            app.figTrajP0.Position = [50 50 925 300];
            app.figTrajP0.Name = 'p0';
            
            t=tiledlayout(app.figTrajP0,1,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axTrajP0 = nexttile(t);
%             app.axTrajP0.XGrid='on';
%             app.axTrajP0.YGrid='on';
            xlabel(app.axTrajP0,'t')
            ylabel(app.axTrajP0,'x')
%             app.axTrajP0 = axes(app.figTrajP0);
            
            app.figTrajP0.Visible = 'on';
        end
        
        function createTrajClickFig(app)
            % Create figTraj and hide until all components are created
            app.figTrajClick = figure('Visible', 'off');
            app.figTrajClick.Position = [1000 50 900 950];
            app.figTrajClick.Name = 'click trajectories';

            app.tilesTrajClick=tiledlayout(app.figTrajClick,app.nClick,1);
            app.tilesTrajClick.TileSpacing = 'compact';
            app.tilesTrajClick.Padding = 'compact';
            for i=1:app.nClick
                app.axTrajClick(i) = nexttile(app.tilesTrajClick);
                app.axTrajClick(i).XTickLabel=[];
%                 app.axTrajClick(i).XAxis.Visible='off';
%                 app.axTrajClick(i).XGrid='on';
%                 app.axTrajClick(i).YGrid='on';
            end
            app.axTrajClick(end).XTickLabelMode='auto';
            xlabel(app.axTrajClick(3),'t')
            ylabel(app.axTrajClick(3),'x')

            app.figTrajClick.Visible = 'on';
        end
    end
    
end