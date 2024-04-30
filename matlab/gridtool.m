classdef gridtool < handle %matlab.mixin.SetGet
    %GRIDTOOL is a GUI for exploring the dependence of solutions of ODEs on
    %parameters, two-parameter planes at a time.  The ODEs are solved in
    %parallel using OpenCL via MEX interface to the c++ clODE library. The
    %interface is modeled after Bard Ermentrout's XPPAUT: ODE files that
    %run in XPP are parsed and converted to the required OpenCL code
    %automatically, and the interface shares some of the same keyboard
    %shortcuts.
    
    %TODO: store the app object as User Data in Control fig?

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
        z %struct with z info: row from gridvars
        dz %increment to change z when using +/-
        XF %place to store the final state of most recent grid simulation
        
        gridTspan
        trajTspan
        
        doContours=false

%         gridSol=struct('x',[],'y',[],'F',[]);
        %mechanism to cache F(:,:,z) values - for rapidly scanning z post compute
        
        feature %map specifying feature functions for display
        fscale=1
        tscale=1

        nClick=3
        trajP0=struct('p0',[],'x0',[],'tspan',[],...
            't',[],'x',[],'dx',[],'aux',[],'xp',nan,'yp',nan,'zp',nan); %struct for trajectory data
        trajClick=struct('p0',[],'x0',[],'tspan',[],...
            't',[],'x',[],'dx',[],'aux',[],'xp',nan,'yp',nan,'zp',nan); %struct for trajectory data
    end
    
    properties (Dependent = true)
        gridx %grid xcoords
        gridy %grid ycoords
        p0 %current parameter vector value (extract from gridvars)
        x0 %current initial condition value
        P0 %base parameter matrix with p0 replicated nPts tim es [nPts, nPar]
        X0 %base ic matrix with x0 replicated nPts times
        nPts
    end
    
    %data needed to display graphs: set-observable, so that whenever these
    %change, triggers update to relevant figures
    properties (SetObservable = true)
        
        gridvars %table of possible grid vars (par+var) with values
        gridtab=table('Size',[2,6], 'VariableTypes',...
            {'cellstr','double','double','double','categorical','double'},...
            'VariableNames',{'name','val','lb','ub','type','ix'},...
            'RowNames',{'x','y'}) %table specifying grid x, y info
    end
    
    %listeners
    properties (Access = private)
        listenerGrid %call updateGridChanges
        listenerGridVars %call updateGridChanges
        listenerFeature
        listenerGridSol
        listenerTrajP0Sol
        listenerTrajClickSol
    end
    
    % Properties that correspond to UI components
    properties (Access = public)
        
        %UI state
        gridclBuilt=false %flag to indicate if OpenCL program is built
        trajclBuilt=false
        gridclInitialized=false %flag to indicate if ready for integration
        trajclInitialized=false
        
        gridIsCurrent=false  %flag to indicate need for new solution
        trajP0SolIsCurrent=false
        trajClickIsCurrent=false
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %grid figure
        figGrid                 matlab.ui.Figure
        axGrid                  matlab.graphics.axis.Axes
        imGrid                  matlab.graphics.primitive.Image
        gridCBar
        markerP0                matlab.graphics.primitive.Line
        markerPquery            matlab.graphics.primitive.Line
        
        gridContours            matlab.graphics.chart.primitive.Contour
        gridColormap
        cLevs = 7

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %figure with p0 trajectory
        figTrajP0                 matlab.ui.Figure
        axyyTrajP0                matlab.graphics.axis.Axes
        lineyyTrajP0              matlab.graphics.chart.primitive.Line 
        ax3DTrajP0                matlab.graphics.axis.Axes
        line3DTrajP0              matlab.graphics.chart.primitive.Line 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %click trajectory figure with n trajectories
        figTrajClick              matlab.ui.Figure
        tilesTrajClick
        axyyTrajClick             matlab.graphics.axis.Axes
        lineyyTrajClick           matlab.graphics.chart.primitive.Line 
        ax3DTrajClick             matlab.graphics.axis.Axes
        line3DTrajClick           matlab.graphics.chart.primitive.Line 
        
        hlinkyyL             matlab.graphics.internal.LinkProp
        hlinkyyR             matlab.graphics.internal.LinkProp
        hlink3d             matlab.graphics.internal.LinkProp
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %control figure
        %TODO: gridlayout instead of manual pixel positions everywhere
        
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
        
        SolverParLabel          matlab.ui.control.Label
        SolverParTable          matlab.ui.control.Table
        gridObserverParTable    matlab.ui.control.Table
        gridObserverParLabel    matlab.ui.control.Label
        
        GridTab                 matlab.ui.container.Tab
        
        gridTable               matlab.ui.control.Table
        gridZDropDownLabel      matlab.ui.control.Label
        gridZDropDown           matlab.ui.control.DropDown
        dzEditFieldLabel        matlab.ui.control.Label
        dzEditField             matlab.ui.control.NumericEditField
        zValueEditField         matlab.ui.control.NumericEditField
        tspanTable              matlab.ui.control.Table
        featureDropDownLabel    matlab.ui.control.Label
        featureDropDown         matlab.ui.control.DropDown
        fscaleEditFieldLabel    matlab.ui.control.Label
        fscaleEditField         matlab.ui.control.NumericEditField
        tscaleEditFieldLabel    matlab.ui.control.Label
        tscaleEditField         matlab.ui.control.NumericEditField
        parLabel                matlab.ui.control.Label
        parTable                matlab.ui.control.Table
        icLabel                 matlab.ui.control.Label
        icTable                 matlab.ui.control.Table
        parDefaultButton        matlab.ui.control.Button
        parSaveButton           matlab.ui.control.Button
        icDefaultButton         matlab.ui.control.Button
        
%         TrajTab                 matlab.ui.container.Tab
        
        trajectoriesLabel       matlab.ui.control.Label
        showTrajP0CheckBox      matlab.ui.control.CheckBox
        showTrajClickCheckBox   matlab.ui.control.CheckBox
        clicksEditFieldLabel    matlab.ui.control.Label
        clicksEditField         matlab.ui.control.NumericEditField
        trajXDropDownLabel      matlab.ui.control.Label
        trajXDropDown           matlab.ui.control.DropDown
        trajYDropDownLabel      matlab.ui.control.Label
        trajYDropDown           matlab.ui.control.DropDown
        trajZDropDownLabel      matlab.ui.control.Label
        trajZDropDown           matlab.ui.control.DropDown
        linkAxesButton          matlab.ui.control.StateButton
        threeDButton            matlab.ui.control.StateButton
        
        %axis in focus?
        
%         HelpTab                 matlab.ui.container.Tab
    end
    
    
    methods (Access = public)
        %constructor
        function app=gridtool(odefile, opts)
            arguments
                odefile=[]
                opts.observer='basicall'
                opts.stepper='dopri5'
                opts.precision='single'
                opts.nGrid=[32,32]
                opts.nClick=3
                opts.tscale=1
                opts.grid_device="GPU"
                opts.traj_device="CPU"
                opts.colormap = [0.9*[1,1,1];turbo]
            end
            
            app.gridColormap=opts.colormap;

            app.createControlFig();
            app.createGridFig();
            app.createTrajP0Fig();
            app.createTrajClickFig();
            app.attachListeners();
                        
            %auto select first GPU for grid, fastest clock for traj
            app.devices=queryOpenCL();
            app.grid_device=find({app.devices(:).type}==opts.grid_device,1,'first');
            app.traj_device=find({app.devices(:).type}==opts.traj_device,1,'first');
%             [~,app.traj_device]=max([app.devices(:).maxClock]);
            
            app.gridDeviceDropDown.Items={app.devices(:).name};
            app.gridDeviceDropDown.ItemsData=1:length(app.devices);
            app.gridDeviceDropDown.Value=app.grid_device;
            app.trajDeviceDropDown.Items={app.devices(:).name};
            app.trajDeviceDropDown.ItemsData=1:length(app.devices);
            app.trajDeviceDropDown.Value=app.traj_device;

            app.nGrid=opts.nGrid;
            app.nClick=opts.nClick;
            app.tscale=opts.tscale;
            app.tscaleEditField.Value=opts.tscale;

            app.loadNewProblem(odefile)

            app.clo_g=clODEfeatures(app.prob, opts.precision, app.grid_device, opts.stepper, opts.observer);
            app.clo_t=clODEtrajectory(app.prob, opts.precision, app.traj_device, opts.stepper);

            app.processNewProblem()

            app.gridObserverDropDown.Items=app.clo_g.observerNames;
            app.gridObserverDropDown.Value=opts.observer;
            app.gridStepperDropDown.Items=app.clo_g.getAvailableSteppers();
            app.gridStepperDropDown.Value=opts.stepper;
            app.trajStepperDropDown.Items=app.clo_t.getAvailableSteppers();
            app.trajStepperDropDown.Value=opts.stepper;
            
            defaultsp=clODE.defaultSolverParams(); %app.clo_g.sp
            spvalues=cell2mat(struct2cell(defaultsp));
            sptable=table('RowNames',fieldnames(defaultsp));
            sptable.grid=spvalues;
            sptable.traj=spvalues;
            app.SolverParTable.Data=sptable;

            %solver opts from ODE file
            % needs check for presence
            app.SolverParTable.Data{'dt',:}=app.prob.opt.dt;
            app.SolverParTable.Data{'dtmax',:}=app.prob.opt.dtmax;
            app.SolverParTable.Data{'max_store',:}=app.prob.opt.maxstor;
            app.SolverParTable.Data{'abstol',:}=app.prob.opt.atoler;
            app.SolverParTable.Data{'reltol',:}=app.prob.opt.toler;
            
            defaultop=clODEfeatures.defaultObserverParams(); %app.clo_g.op
            optable=table('RowNames',fieldnames(defaultop));
            optable.value=cell2mat(struct2cell(defaultop));
            app.gridObserverParTable.Data=optable;
            app.gridObserverParTable.RowName=fieldnames(defaultop);

            %             if nargout==0
            %                 clear app
            %             end
            
        end
        
        % Copy current parameter set into an XPP file (requires
        % ChangeXPPodeFile and packagePars4XPP)
        function savePars(app,~,~)
            
            [fileToWrite,path]=uiputfile('*.ode','Save parameters to ODE file', app.odefile);
            
            fileToWrite=fullfile(path,fileToWrite);
            
            if ~exist(fileToWrite,'file')
                copyfile(app.odefile,fileToWrite)
            end
            
            ChangeXPPodeFile(fileToWrite,package4XPP(app.prob.parNames, app.p0, app.prob.varNames, app.x0));
        end
        
        function loadNewProblem(app,odefile)
            [~,app.prob]=ode2cl(odefile);
            app.odefile=app.prob.name+".ode";
        end
        
        function processNewProblem(app)
            
            if app.prob.nPar+app.prob.nVar<2
                error('Not enough parameters and variables for a 2D grid!')
            end
            
            app.currentFileLabel.Text=app.prob.name+".ode";
            
            if isempty(app.clo_g)
                
            else %to avoid reverting to default sp and op:
                app.clo_g.setProblemInfo(app.prob);
                app.clo_t.setProblemInfo(app.prob);
                
                %remove grid.name to allow new grid to be set below
                app.listenerGrid.Enabled=0;
                app.gridtab.name={'';''};
                app.listenerGrid.Enabled=1;
                app.clo_g.F=[]; %to allow new nVar
                app.clo_g.Xf=[]; %to allow new nVar
            end        
            
            app.fscale=1;
            app.XF=[]; %to allow new nVar
                
            %tspan
            t0=app.prob.opt.t0;
            tf=app.prob.opt.total;
            app.clo_g.tspan=[t0,tf];
            app.clo_t.tspan=[t0,tf];
            app.gridTspan=[t0,tf];
            app.trajTspan=[t0,tf];
            
            tspan=table;
            tspan.t0=[t0;t0];
            tspan.tf=[tf;tf];
%             tspan.trans=[tf;tf];
            app.tspanTable.Data=tspan;
            
            %valid grid variables
            icNames=strcat(app.prob.varNames(:),'0');
            newGridvars=table;
            newGridvars.name=[app.prob.parNames(:);icNames(:)];
            newGridvars.val=[app.prob.p0(:);app.prob.x0(:)];
            newGridvars.lb=[[app.prob.par.lb]';[app.prob.var.lb]'];
            newGridvars.ub=[[app.prob.par.ub]';[app.prob.var.ub]'];
            newGridvars.type=categorical([ones(app.prob.nPar,1);2*ones(app.prob.nVar,1)],[1,2],{'par','ic'});
            newGridvars.ix=[(1:app.prob.nPar)';(1:app.prob.nVar)'];
            newGridvars.Properties.RowNames=newGridvars.name;
            
            %autobounds - should make this a button
            const_bounds=newGridvars.lb==newGridvars.ub;
%             const_bounds=const_bounds & newGridvars.type=="par"; %don't do ic?
            const_vals=newGridvars.val(const_bounds);
            newGridvars.lb(const_bounds)=const_vals-abs(const_vals)*0.5;
            newGridvars.ub(const_bounds)=const_vals+abs(const_vals)*0.5;
            
            %set grid to first two names
            app.gridtab.name=newGridvars.name(1:2);
            
            % specify Z for +/- incrementing
            app.z=newGridvars(3,:);
            app.dz=(app.z.ub-app.z.lb)/10;
            app.gridZDropDown.Items=newGridvars.name;
            app.gridZDropDown.Value=app.z.name;
            app.dzEditField.Value=app.dz;
            app.zValueEditField.Value=app.z.val;
            
            app.gridvars=newGridvars; %triggers postSet listener

            %set up selectable gridTable name
            app.gridTable.ColumnFormat={app.gridvars.name',[],[],[],[]};
            
            % trajectory variables
            app.trajvars=[app.prob.varNames(:);app.prob.auxNames(:);...
                strcat('d',app.prob.varNames(:),'/dt')];
            
            %trajectory plot variables
            app.trajXDropDown.Items=[{'t'};app.trajvars];
            app.trajYDropDown.Items=app.trajvars;
            app.trajZDropDown.Items=[{'none'};app.trajvars];
            
            yyaxis(app.axyyTrajP0,'left');
            xlabel(app.axyyTrajP0,app.trajXDropDown.Value);
            ylabel(app.axyyTrajP0,app.trajYDropDown.Value);
            yyaxis(app.axyyTrajClick(end),'left');
            xlabel(app.axyyTrajClick(end),app.trajXDropDown.Value);
            ylabel(app.axyyTrajClick(end),app.trajYDropDown.Value);
           
            app.buildCL('grid');
            app.buildCL('traj');
            
            figure(app.figControl)
        end
    
        
        function buildCL(app,which_clo)
            switch which_clo
                case 'grid'
                    app.clo_g.buildCL();
                    app.updateFeatureNames();
                    app.clo_g.Xf=app.X0;
                    app.makeGridData();
                    app.clo_g.initialize();
                    app.gridCLIsBuiltButton.Value = 1;
                    app.gridCLIsBuiltButton.Text = 'CL ready';
                    app.gridCLIsBuiltButton.BackgroundColor=[0.7,1,0.7];
                case 'traj'
                    app.clo_t.buildCL();
                    app.clo_t.P=app.p0;
                    app.clo_t.X0=app.x0;
                    app.clo_t.Xf=app.x0;
                    app.makeTrajData(); %p0: no clickCoords
                    app.clo_t.initialize();
                    app.trajCLIsBuiltButton.Value = 1;
                    app.trajCLIsBuiltButton.Text = 'CL ready';
                    app.trajCLIsBuiltButton.BackgroundColor=[0.7,1,0.7];
            end
        end
        
        function updateFeatureNames(app)
            app.clo_g.getNFeatures(); %update features (nVar may affect)
            app.clo_g.featureNames();
            fNames=app.clo_g.featureNames();
            ampvars=fNames(startsWith(fNames,'max')); %same as trajvars
            ampvars=split(ampvars);
            ampvars=ampvars(:,2);
            f=containers.Map;
            fNamesPlus=string();
            for i=1:length(ampvars)
                var=ampvars{i};
                f(var+" max")=@(F)F(:,fNames=="max "+var);
                fNamesPlus(end+1,1)=var+" max";
                f(var+" min")=@(F)F(:,fNames=="min "+var);
                fNamesPlus(end+1,1)=var+" min";
                if any(fNames=="mean "+var)
                    f(var+" mean")=@(F)F(:,fNames=="mean "+var);
                    fNamesPlus(end+1,1)=var+" mean";
                end
                f(var+" range")=@(F)F(:,fNames=="max "+var)-F(:,fNames=="min "+var);
                fNamesPlus(end+1,1)=var+" range";
            end
            extraF=fNames(contains(fNames,'count'));
            for i=1:length(extraF)
                thisF=string(extraF{i});
                f(thisF)=@(F)F(:,fNames==thisF);
            end
            fNamesPlus(1)=[];
            fNamesPlus=[fNamesPlus;fNames(contains(fNames,'count'))];
            app.feature=f;
            app.featureDropDown.Items=fNamesPlus;
            app.featureDropDown.Value=fNamesPlus(1);
            title(app.gridCBar,app.featureDropDown.Value);
        end
            
    end



    % Callbacks that handle component events
    methods (Access = private)
        
        % Close request function: figControl
        function figControlCloseRequest(app, src, event)
            if isvalid(app.figGrid), close(app.figGrid); end
            if isvalid(app.figTrajP0), close(app.figTrajP0); end
            if isvalid(app.figTrajClick), close(app.figTrajClick); end
            closereq
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
            app.loadNewProblem(app.odefile)
            app.processNewProblem()
        end
        
        % Value changed function: gridDeviceDropDown
        function gridDeviceDropDownValueChanged(app, src, event)
            value = app.gridDeviceDropDown.Value;
            app.clo_g.setOpenCL(value);
            app.setBuildNeeded('grid');
            app.gridclInitialized=false;
        end
        
        % Value changed function: trajDeviceDropDown
        function trajDeviceDropDownValueChanged(app, src, event)
            value = app.trajDeviceDropDown.Value;
            app.clo_t.setOpenCL(value);
            app.setBuildNeeded('traj');
            app.trajclInitialized=false;
        end
        
        % Value changed function: gridPrecisionCheckBox
        function gridPrecisionCheckBoxValueChanged(app, src, event)
            useDouble = app.gridPrecisionCheckBox.Value;
            precision='single';
            if useDouble, precision='double'; end
            app.clo_g.setPrecision(precision);
            app.setBuildNeeded('grid');
            app.gridclInitialized=false;
        end
        
        % Value changed function: trajPrecisionCheckBox
        function trajPrecisionCheckBoxValueChanged(app, src, event)
            useDouble = app.trajPrecisionCheckBox.Value;
            precision='single';
            if useDouble, precision='double'; end
            app.clo_t.setPrecision(precision);
            app.setBuildNeeded('traj');
            app.trajclInitialized=false;
        end
        
        % Value changed function: gridStepperDropDown
        function gridStepperDropDownValueChanged(app, src, event)
            stepper = app.gridStepperDropDown.Value;
            app.clo_g.setStepper(stepper);
            app.setBuildNeeded('grid');
            app.gridclInitialized=false;
        end
        
        % Value changed function: trajStepperDropDown
        function trajStepperDropDownValueChanged(app, src, event)
            stepper = app.trajStepperDropDown.Value;
            app.clo_t.setStepper(stepper);
            app.setBuildNeeded('traj');
            app.trajclInitialized=false;
        end
        
        % Value changed function: gridObserverDropDown
        function gridObserverDropDownValueChanged(app, src, event)
            observer = app.gridObserverDropDown.Value;
            app.clo_g.setObserver(observer);
            app.clo_g.getNFeatures(); %update features (nVar may affect)
            app.featureDropDown.Items=app.clo_g.featureNames();
            app.setBuildNeeded('grid');
            app.gridclInitialized=false;
        end
        
        function augmentFeatures(app)
            basefeatures=app.clo_g.featureNames();
            
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
        
        function setBuildNeeded(app,which_clo)
            %listen for changes in clo.clBuilt? Two ifs: if
            %~clo_g.clBuilt, if ~clo_t.clBuilt?
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
        
        
        % Cell edit callback: SolverParTable
        function SolverParTableCellEdit(app, src, event)
            indices = event.Indices;
            newData = event.NewData;
            thisfield=src.DisplayData.Properties.RowNames{indices(1)};
            if isnumeric(newData) && newData>0
                if indices(2)==1
                    sp=app.clo_g.sp;
                    sp.(thisfield)=newData;
                    app.clo_g.setSolverPars(sp); %do the device transfer (could clean this up if clODE had set.sp)
                elseif indices(2)==2
                    sp=app.clo_t.sp;
                    sp.(thisfield)=newData;
                    app.clo_t.setSolverPars(sp);
                end
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
        
        
        % Cell edit callback: gridTable
        function gridTableCellEdit(app, src, event)
            % gridTable val/lb/ub must be synced to parTable and icTable
            indices = event.Indices;
            newData = event.NewData;
            vix=indices(1);
            prop=src.ColumnName{indices(2)};
            switch prop
                case 'name'
                    app.gridtab(vix,:)=app.gridvars(newData,:);
                    newGridDisplayData=app.gridtab(vix,{'name','val','lb','ub'});
                    newGridDisplayData.N=app.nGrid(vix);
                    app.gridTable.Data(vix,:)=table2cell(newGridDisplayData);
                    app.updateGridPlot('nosol');
                    app.makeGridData(); %prep for next integration
                    
                case {'val', 'lb', 'ub'}
                    app.gridvars{app.gridtab.name(vix),prop}=newData; %sync gridvars table
                    
                case 'N'
                    app.nGrid(indices(1))=newData;
                    app.updateGridPlot('nosol');
                    app.makeGridData(); %prep for next integration
            end
        end
        
        % Value changed function: gridZDropDown
        function gridZDropDownValueChanged(app, src, event)
            newZ = app.gridZDropDown.Value;
            app.z=app.gridvars(newZ,:);
            app.dz=(app.z.ub-app.z.lb)/10;
            app.dzEditField.Value=app.dz;
            app.zValueEditField.Value=app.z.val;
        end
        
        % Value changed function: dzEditField
        function dzEditFieldValueChanged(app, src, event)
            newDZ = app.dzEditField.Value;
            if isnumeric(newDZ) && isreal(newDZ)
                app.dz=newDZ;
            end
        end
        
        % Value changed function: zValueEditField
        function zValueEditFieldValueChanged(app, src, event)
            newZval = app.zValueEditField.Value;
            if isnumeric(newZval) && isreal(newZval)
                app.gridvars{app.z.name,'val'}=newZval;
            end
        end
        
        
        % Cell edit callback: tspanTable
        function tspanTableCellEdit(app, src, event)
            indices = event.Indices; %r=1:grid, r=2:traj, c=1:t0, c2:tf, c3:trans
            newData = event.NewData;
            if indices(1)==1
                app.gridTspan(indices(2))=newData;
                app.clo_g.settspan(app.gridTspan);
            elseif indices(1)==2
                app.trajTspan(indices(2))=newData;
                app.clo_t.settspan(app.trajTspan);
            end
        end
        
        
        % Value changed function: featureDropDown
        function featureDropDownValueChanged(app, src, event)
            newFeature = app.featureDropDown.Value;
            title(app.gridCBar,newFeature,'Interpreter','none')
            app.updateGridPlot('feature');
        end
        
        % Value changed function: fscaleEditField
        function fscaleEditFieldValueChanged(app, src, event)
            app.fscale = app.fscaleEditField.Value;
            app.updateGridPlot('feature');
        end
        
        % Value changed function: tscaleEditField
        function tscaleEditFieldValueChanged(app, src, event)
            %change t scale for trajectories (only?)
            %TODO: change should really be sent to feature detector, so
            %time-dep features are also changed?
            app.tscale = app.tscaleEditField.Value;
%             app.trajP0.t=app.trajP0.t*app.tscale;
            app.updateTrajPlotData('p0');
            app.updateTrajPlotData('click');
            app.updateTrajPlot('p0');
            app.updateTrajPlot('click');
        end
        
        
        % Cell edit callback: parTable
        function parTableCellEdit(app, src, event)
            %update clo_g.P. Link any lb/ub changes to gridTable.
            indices = event.Indices;
            newData = event.NewData;
            var=src.Data.name(indices(1));
            prop=src.ColumnName{indices(2)};
            app.gridvars{var, prop}=newData; %sync gridvars table
            if prop=="val"
                app.makeGridData(); %non-grid param value was updated
            end
        end
        
        % Cell edit callback: icTable
        function icTableCellEdit(app, src, event)
            %update clo_g.X0. Link any lb/ub changes to gridTable.
            indices = event.Indices;
            newData = event.NewData;
            var=src.Data.name(indices(1));
            prop=src.ColumnName{indices(2)};
            app.gridvars{var, prop}=newData; %sync gridvars table
        end
        
        % Button pushed function: parSaveButton
        function parSaveButtonPushed(app, src, event)
            app.savePars();
        end
        
        % Button pushed function: parDefaultButton
        function parDefaultButtonPushed(app, src, event)
            %get the original default values for par.val from app.prob
            newGridvars=table;
            newGridvars.val=app.prob.p0(:);
            app.gridvars(app.gridvars.type=="par",{'val'})=newGridvars;
            %get the original default values for par.val/lb/ub from app.prob
%             newGridvars.lb=[app.prob.par.lb]';
%             newGridvars.ub=[app.prob.par.ub]';
%             const_bounds=newGridvars.lb==newGridvars.ub;
%             const_vals=newGridvars.val(const_bounds);
%             newGridvars.lb(const_bounds)=const_vals-abs(const_vals)*0.05;
%             newGridvars.ub(const_bounds)=const_vals+abs(const_vals)*0.05;
%             app.gridvars(app.gridvars.type=="par",{'val','lb','ub'})=newGridvars;
            app.makeGridData(); %in case grid wasn't changed, but other pars were
        end
        
        % Button pushed function: icDefaultButton
        function icDefaultButtonPushed(app, src, event)
            %get the original default values for var.val/lb/ub from app.prob
            newGridvars=table;
            newGridvars.val=app.prob.x0(:);
            app.gridvars(app.gridvars.type=="ic",{'val'})=newGridvars;
%             newGridvars.lb=[app.prob.var.lb]';
%             newGridvars.ub=[app.prob.var.ub]';
%             const_bounds=newGridvars.lb==newGridvars.ub;
%             const_vals=newGridvars.val(const_bounds);
%             newGridvars.lb(const_bounds)=const_vals-abs(const_vals)*0.05;
%             newGridvars.ub(const_bounds)=const_vals+abs(const_vals)*0.05;
%             app.gridvars(app.gridvars.type=="ic",{'val','lb','ub'})=newGridvars;
            app.makeGridData(); %in case grid wasn't changed, but other pars were
        end
        
        %CheckBox value changed function: showTrajP0CheckBox
        function showTrajP0CheckBoxValueChanged(app, src, event)
            if app.showTrajP0CheckBox.Value
                app.markerP0.Visible='on';
                app.figTrajP0.Visible='on';
            else
                app.markerP0.Visible='off';
                app.figTrajP0.Visible='off';
            end
        end
        
        
        %CheckBox value changed function: showTrajClickCheckBox
        function showTrajClickCheckBoxValueChanged(app, src, event)
            if app.showTrajClickCheckBox.Value
                for i=1:app.nClick
                    app.markerPquery(i).Visible='on';
                end
                app.figTrajClick.Visible='on';
            else
                for i=1:app.nClick
                    app.markerPquery(i).Visible='off';
                end
                app.figTrajClick.Visible='off';
            end
        end
        
        % Value changed function: trajXDropDown
        function trajXDropDownValueChanged(app, src, event)
            newX = app.trajXDropDown.Value;
            app.updateTrajPlotData('p0');
            app.updateTrajPlotData('click');
            app.updateTrajPlot('p0');
            app.updateTrajPlot('click');
        end
        
        % Value changed function: trajYDropDown
        function trajYDropDownValueChanged(app, src, event)
            newY = app.trajYDropDown.Value;
            app.updateTrajPlotData('p0');
            app.updateTrajPlotData('click');
            app.updateTrajPlot('p0');
            app.updateTrajPlot('click');
        end
        
        % Value changed function: trajZDropDown
        function trajZDropDownValueChanged(app, src, event)
            app.updateTrajPlotData('p0');
            app.updateTrajPlotData('click');
            app.updateTrajPlot('p0');
            app.updateTrajPlot('click');
        end
        
        % Value changed function: clicksEditField
        function clicksEditFieldValueChanged(app, src, event)
            newNClick = app.clicksEditField.Value;
            app.nClick=newNClick;
            app.createTrajClickFig(); %don't make new Traj Data: needs new clicks
        end
        
        % Value changed function: linkAxesButton
        function linkAxesButtonValueChanged(app, src, event)
            if app.linkAxesButton.Value
                app.setLinkedLims();
%                 axyyL=app.axyyTrajP0.YAxis(1);
%                 axyyR=app.axyyTrajP0.YAxis(2);
%                 for i=1:app.nClick
%                     axyyL=[axyyL;app.axyyTrajClick(i).YAxis(1)];
%                     axyyR=[axyyR;app.axyyTrajClick(i).YAxis(2)];
%                 end
%                 app.hlinkyyL=linkprop(axyyL,'Limits');
%                 app.hlinkyyR=linkprop(axyyR,'Limits');
%                 
%                 ax3d=app.ax3DTrajP0;
%                 for i=1:app.nClick
%                     ax3d=[ax3d;app.ax3DTrajClick(i)];
%                 end
%                 app.hlink3d=linkprop(ax3d,{'XLim','YLim','ZLim'});
                
            else
%                 app.hlinkyyL.removeprop('Limits')
%                 app.hlinkyyR.removeprop('Limits')
%                 app.hlink3d.removeprop('XLim')
%                 app.hlink3d.removeprop('YLim')
%                 app.hlink3d.removeprop('ZLim')
                
                %yyaxes
                app.axyyTrajP0.XLimMode='auto';
                app.axyyTrajP0.YAxis(1).LimitsMode='auto';
                app.axyyTrajP0.YAxis(2).LimitsMode='auto';
                for i=1:app.nClick
                    app.axyyTrajClick(i).XLimMode='auto';
                    app.axyyTrajClick(i).YAxis(1).LimitsMode='auto';
                    app.axyyTrajClick(i).YAxis(2).LimitsMode='auto';
                end
                %3daxes
                app.ax3DTrajP0.XLimMode='auto';
                app.ax3DTrajP0.YLimMode='auto';
                app.ax3DTrajP0.ZLimMode='auto';
                for i=1:app.nClick
                    app.ax3DTrajClick(i).XLimMode='auto';
                    app.ax3DTrajClick(i).YLimMode='auto';
                    app.ax3DTrajClick(i).ZLimMode='auto';
                end
            end
        end
        
        function setLinkedLims(app)
            
            xL=[app.lineyyTrajP0(1).XData,app.lineyyTrajClick(:,1).XData];
            xLimL=[min(xL),max(xL)]; 
            if diff(xLimL)==0
                xLimL=xLimL(1)+[-1,1]*0.1*xLimL(1);
            else
                xLimL=xLimL+[-1,1]*diff(xLimL)*0.01;
            end

            yL=[app.lineyyTrajP0(1).YData,app.lineyyTrajClick(:,1).YData];
            yLimL=[min(yL),max(yL)]; 
            if diff(yLimL)==0
                yLimL=yLimL(1)+[-1,1]*0.1*yLimL(1);
            else
                yLimL=yLimL+[-1,1]*diff(yLimL)*0.01;
            end

            if ~isnan(app.trajClick(1).zp)
                yR=[app.lineyyTrajP0(2).YData,app.lineyyTrajClick(:,2).YData];
                yLimR=[min(yR),max(yR)]; 
                if diff(yLimR)==0
                    yLimR=yLimR(1)+[-1,1]*0.1*yLimR(1);
                else
                    yLimR=yLimR+[-1,1]*diff(yLimR)*0.01;
                end
            end

            app.axyyTrajP0.XLim=xLimL;
            app.axyyTrajP0.YAxis(1).Limits=yLimL;
            if ~isnan(app.trajClick(1).zp)
                app.axyyTrajP0.YAxis(2).Limits=yLimR;
            end
            for i=1:app.nClick
                app.axyyTrajClick(i).XLim=xLimL;
                app.axyyTrajClick(i).YAxis(1).Limits=yLimL;
                if ~isnan(app.trajClick(1).zp)
                    app.axyyTrajClick(i).YAxis(2).Limits=yLimR;
                end
            end

            %3daxes
            if ~isnan(app.trajClick(1).zp)
                xx=[app.line3DTrajP0.XData,app.line3DTrajClick(:).XData];
                yy=[app.line3DTrajP0.YData,app.line3DTrajClick(:).YData];
                zz=[app.line3DTrajP0.ZData,app.line3DTrajClick(:).ZData];

                xLim=[min(xx),max(xx)];
                if diff(xLim)==0
                    xLim=xLim(1)+[-1,1]*0.1*xLim(1);
                else
                    xLim=xLim+[-1,1]*diff(xLim)*0.01;
                end
                yLim=[min(yy),max(yy)];
                if diff(yLim)==0
                    yLim=yLim(1)+[-1,1]*0.1*yLim(1);
                else
                    yLim=yLim+[-1,1]*diff(yLim)*0.01;
                end
                zLim=[min(zz),max(zz)];
                if diff(zLim)==0
                    zLim=zLim(1)+[-1,1]*0.1*zLim(1);
                else
                    zLim=zLim+[-1,1]*diff(zLim)*0.01;
                end

                app.ax3DTrajP0.XLim=xLim;
                app.ax3DTrajP0.YLim=yLim;
                app.ax3DTrajP0.ZLim=zLim;
                for i=1:app.nClick
                    app.ax3DTrajClick(i).XLim=xLim;
                    app.ax3DTrajClick(i).YLim=yLim;
                    app.ax3DTrajClick(i).ZLim=zLim;
                end
            end
        end
        
        % Value changed function: threeDButton
        function threeDButtonValueChanged(app, src, event)
            do3D = app.threeDButton.Value;
            if do3D
                app.trajYDropDownLabel.Text = 'y';
                app.trajZDropDownLabel.Text = 'z';
            else
                app.trajYDropDownLabel.Text = 'y1';
                app.trajZDropDownLabel.Text = 'y2';
            end
            app.updateTrajPlotData('p0');
            app.updateTrajPlotData('click');
            app.updateTrajPlot('p0');
            app.updateTrajPlot('click');
        end
    end
    
    
    %get/set
    methods
        function N=get.nPts(app)
            N=prod(app.nGrid);
        end
        
        function p0=get.p0(app)
            p0=app.gridvars.val(app.gridvars.type=="par")';
        end
        
        function x0=get.x0(app)
            x0=app.gridvars.val(app.gridvars.type=="ic")';
        end
        
        function P0=get.P0(app)
            P0=repmat(app.gridvars.val(app.gridvars.type=="par")',app.nPts,1);
        end
        
        function X0=get.X0(app)
            X0=repmat(app.gridvars.val(app.gridvars.type=="ic")',app.nPts,1);
        end
        
        function gridx=get.gridx(app)
            gridx=linspace(app.gridtab.lb(1),app.gridtab.ub(1),app.nGrid(1));
        end
        
        function gridy=get.gridy(app)
            gridy=linspace(app.gridtab.lb(2),app.gridtab.ub(2),app.nGrid(2));
        end
    end
    
    %listeners
    methods
        function attachListeners(app)
%             app.listenerGrid=addlistener(app,'grid','PostSet',@gridtool.gridUpdate);
            app.listenerGridVars=addlistener(app,'gridvars','PostSet',@gridtool.gridvarUpdate);
        end
    end
    
    
%TODO: just use gridTable.Data everywhere. no grid update listener???
    methods (Static = true)
%         function gridUpdate(~,eventData)
%             app = eventData.AffectedObject;
%             newGridDisplayData=app.gridtab(:,{'name','val','lb','ub'});
%             newGridDisplayData.N=app.nGrid(:);
%             app.gridTable.Data=table2cell(newGridDisplayData);
%             app.updateGridPlot('nosol');
%             app.makeGridData(); %prep for next integration
%         end
        
        function gridvarUpdate(~,eventData)
            app = eventData.AffectedObject;
            %update grid
%             if isempty(app.gridtab.name{1})
%                 vars=app.gridvars.name(1:2);
%             else
                vars=app.gridtab.name;
%             end
            newGrid=app.gridvars(vars,:);
            newGrid.Row={'x';'y'};
            
            nameChange=~strcmp(newGrid{:,'name'},app.gridtab{:,'name'});
            boundChange=newGrid{:,{'lb','ub'}}~=app.gridtab{:,{'lb','ub'}};
            if any(nameChange(:))||any(boundChange(:)) %grid solution is invalidated
                app.gridtab=newGrid; 
                newGridDisplayData=app.gridtab(:,{'name','val','lb','ub'});
                newGridDisplayData.N=app.nGrid(:);
                app.gridTable.Data=table2cell(newGridDisplayData);
                app.updateGridPlot('nosol');
                app.makeGridData(); %prep for next integration
            end
            
            p0change=newGrid{:,'val'}~=app.gridtab{:,'val'};
            if any(p0change(:))
                app.gridtab{:,'val'}=newGrid{:,'val'};
                app.gridTable.Data(:,2)=num2cell(newGrid{:,'val'});
                app.markerP0.XData=app.gridtab{'x','val'};
                app.markerP0.YData=app.gridtab{'y','val'};
            end
            
            newZ=app.gridvars(app.z.name,:); %zname doesn't change, but val/lb/ub might
            if ~isequal(newZ, app.z)
                app.z=newZ;
%                 app.dz=(app.z.ub-app.z.lb)/10;
%                 app.dzEditField.Value=app.dz; 
                app.zValueEditField.Value=app.z.val;
            end
            
            %update parameter and ic table
            app.parTable.Data=app.gridvars(app.gridvars.type=="par",{'name','val','lb','ub'});
            app.icTable.Data=app.gridvars(app.gridvars.type=="ic",{'name','val','lb','ub'});
            
        end
    end

    % internal functions (clODE control & ui interaction)
    methods (Access = private)
        
        function commonKeyPress(app,src,event)
            disp(event.Key)
            switch(event.Key)
                case 'f1'
                    figure(app.figControl)
                case 'f2'
                    figure(app.figGrid)
                case 'f3'
                    figure(app.figTrajClick)
                case 'f4'
                    figure(app.figTrajP0)
                otherwise
                    if src==app.figGrid
                        app.gridKeyPress(src,event);
                    elseif src~=app.figControl
                        app.trajKeyPress(src,event);
                    end
            end
        end

        function gridKeyPress(app,src,event)
            %             app.clo_g.initialize();
%             disp(event.Key)
            switch(event.Key)
                case 'c' %'continue' - features without initialization
                    app.integrateGrid('continue')
                    
                case 'g' %'go' - features with initialization
                    app.integrateGrid('go')
                    
                case 'p' %'p0' - run the p0 trajectory
                    app.makeTrajData()
                    app.integrateTraj('go','p0')
                    
                case '1' %'random'
                    app.integrateGrid('point')

                case 'r' %'random'
                    app.integrateGrid('random')
                    
                case 't' %'transient'
                    app.integrateGrid('transient')

                case 'z' %zoom grid
                    tmp=get(app.figGrid,'KeyPressFcn'); %
                    set(app.figGrid,'KeyPressFcn',[]); %temporarily turn off clicktrajectory
                    
                    k=waitforbuttonpress;
                    if k==0
                        point1 = app.figGrid.CurrentPoint;    % location when mouse clicked
                        rbbox;                 % rubberbox
                        point2 = app.figGrid.CurrentPoint;    % location when mouse released
                        point1 = point1(1,1:2);            % extract x and y
                        point2 = point2(1,1:2);
                        lb = min(point1,point2);           % calculate bounds
                        ub = max(point1,point2);           % calculate bounds
            
                        if all(ub-lb>0)
                            app.gridvars{app.gridtab.name(1:2),'lb'}=lb;
                            app.gridvars{app.gridtab.name(1:2),'ub'}=ub;
                            app.makeGridData();
                        else
                            disp('Try selecting new bounds more slowly...')
                        end
                    else
                        disp('Window Zoom-in canceled by keypress')
                    end
                    
                    set(app.figGrid,'KeyPressFcn',tmp);
                    
                case 'add' % increment z by +dz, transient
                    newZ=min(app.z.val+app.dz, app.z.ub);
                    if newZ~=app.z.val
                        app.gridvars{app.z.name,'val'}=newZ;
                        app.makeGridData(); %prep for next integration
%                         app.integrateGrid('transient')
%                         app.integrateGrid('random')
%                         app.integrateGrid('point')
                        app.integrateGrid('pointgo')
                    end
                    
                case 'subtract' % increment z by -dz, transient
                    newZ=max(app.z.val-app.dz, app.z.lb);
                    if newZ~=app.z.val
                        app.gridvars{app.z.name,'val'}=newZ;
                        app.makeGridData(); %prep for next integration
%                         app.integrateGrid('transient')
%                         app.integrateGrid('random')
%                         app.integrateGrid('point')
                        app.integrateGrid('pointgo')
                    end
            end
        end
        
        function makeGridData(app)
            newP=app.P0;
            if isempty(app.XF)||numel(app.XF)~=numel(app.X0)
                app.XF=app.X0;
            end
            newX0=app.XF; %defaults to using most recent final state
            [X,Y]=ndgrid(app.gridx, app.gridy);
            
            if app.gridtab{'x','type'}=="par"
                newP(:,app.gridtab{'x','ix'})=X(:);
            else
                newX0(:,app.gridtab{'x','ix'})=X(:);
            end
            if app.gridtab{'y','type'}=="par"
                newP(:,app.gridtab{'y','ix'})=Y(:);
            else
                newX0(:,app.gridtab{'y','ix'})=Y(:);
            end
            
            app.clo_g.setProblemData(newX0, newP);
            figure(app.figGrid)
        end
        
        
        function integrateGrid(app, action)
            %assume all device variables are populated
            tic
            switch action
                case 'continue'
                    app.clo_g.shiftX0(); %device X0<-XF
                    app.clo_g.features(0);
                    app.clo_g.getF();
                    
                case 'go'
                    app.clo_g.shiftX0();
                    app.clo_g.features(1);
                    app.clo_g.getF();
                    
%                 case 'shift'
                case 'point'
                    if app.gridtab{'x','type'}=="par" && app.gridtab{'y','type'}=="par"
                        newX0=repmat(app.gridvars.val(app.gridvars.type=="ic")',app.nPts,1);
                        app.clo_g.setX0(newX0);
                    end
                    app.clo_g.transient();

                case 'pointgo'
                    if app.gridtab{'x','type'}=="par" && app.gridtab{'y','type'}=="par"
                        newX0=repmat(app.gridvars.val(app.gridvars.type=="ic")',app.nPts,1);
                        app.clo_g.setX0(newX0);
                    end
%                     app.clo_g.shiftX0();
                    app.clo_g.transient();
                    app.clo_g.shiftX0();
                    app.clo_g.features(1);
                    app.clo_g.getF();
                    
                case 'random'
                    x0lb=app.gridvars.lb(app.gridvars.type=="ic")';
                    x0ub=app.gridvars.ub(app.gridvars.type=="ic")';
                    newX0=x0lb+rand(app.clo_g.nPts,length(x0lb)).*(x0ub-x0lb);
                    app.clo_g.setX0(newX0);
                    app.clo_g.transient();
                    
                case 'transient'
                    app.clo_g.shiftX0();
                    app.clo_g.transient();
            end
            app.XF=app.clo_g.getXf();
            toc
            
            % call for plot update (use listener to SOL properties?)
            switch action
                case {'continue','go','pointgo'}
                    app.updateGridPlot('feature')
                case {'random','transient','point'}
                    app.updateGridPlot('transient')
            end
        end
        
        function updateGridPlot(app,type)
            if ~isvalid(app.figGrid)
                app.createGridFig();
                app.imGrid.XData=app.gridx;
                app.imGrid.YData=app.gridy;
                app.markerP0.XData=app.gridtab{'x','val'};
                app.markerP0.YData=app.gridtab{'y','val'};
                xlabel(app.axGrid,app.gridtab{'x','name'});
                ylabel(app.axGrid,app.gridtab{'y','name'});
                title(app.gridCBar,app.featureDropDown.Value,'Interpreter','none')
            end
            app.gridContours.Visible='off';
            switch type
                case 'transient'
                    %handle XF plotting selection
                    ix=app.clo_g.op.fVarIx;
                    C=reshape(app.XF(:,ix),app.nGrid)';
                    app.imGrid.CData=C;
                    title(app.gridCBar,app.prob.varNames(ix),'Interpreter','none')
                    
                case 'feature'
                    %handle F plotting seletion
                    if isempty(app.clo_g.F), return; end
                    fun=app.feature(app.featureDropDown.Value);
                    C=reshape(fun(app.clo_g.F)/app.fscale,app.nGrid)';
                    C(C<-1e10|C>1e10)=nan; %hack to prevent display of bad values
                    app.imGrid.CData=C;
                    title(app.gridCBar,app.featureDropDown.Value,'Interpreter','none')

                    if app.doContours && ~sum(C(:))==0
%                         levels = app.gridCBar.Ticks;
                        levels = linspace(round(min(C(:)),1,"significant"), round(max(C(:)),1,"significant"), app.cLevs);
                        app.gridContours.XData = app.gridx;
                        app.gridContours.YData = app.gridy;
                        app.gridContours.ZData = C;
                        app.gridContours.ShowText="on";
                        app.gridContours.LevelList=levels;
%                         hold on
%                         [~,app.gridContours] = contour(app.gridx, app.gridy, C,levels,ShowText="on");
%                         hold off
                        app.gridContours.Visible='on';
                    end
                        
                    
                case 'nosol'  %setup plotting elements when grid is changed
                    app.imGrid.XData=app.gridx;
                    app.imGrid.YData=app.gridy;
                    app.imGrid.CData=zeros(app.nGrid(1),app.nGrid(2));
                    app.markerP0.XData=app.gridtab{'x','val'};
                    app.markerP0.YData=app.gridtab{'y','val'};
                    for i=1:app.nClick
                        app.markerPquery(i).Visible='off';
                    end
                    xlabel(app.axGrid,app.gridtab{'x','name'});
                    ylabel(app.axGrid,app.gridtab{'y','name'});
                    title(app.gridCBar,app.featureDropDown.Value,'Interpreter','none')
            end
            
            if app.imGrid.Visible=="off"
                app.imGrid.Visible='on';
            end
        end
        
        
        function trajKeyPress(app,src,event)
%             disp(event.Key)
            if src==app.figTrajP0
                which_traj='p0';
            elseif src==app.figTrajClick
                which_traj='click';
            end
            
            switch(event.Key)
                case 'c' %'continue' - append+shift
                    app.integrateTraj('continue',which_traj)
                    
                case 'g' %'go' - start here
                    app.integrateTraj('go',which_traj)
                    
                case 'l' %'go' - start here
                    app.integrateTraj('last',which_traj)
                    
                case 'r' %'random'
                    app.integrateTraj('random',which_traj)
                    
                case 's' %'shift'
                    app.integrateTraj('shift',which_traj)
                    
%                 case 't' %'transient'
%                     app.integrateTraj('transient',which_traj)
            end
        end
        
        %clicking on the p0 marker lets the user select new p0 coords
        function clickP0(app,src,event)
            disp('clickP0')
            [px,py]=ginput(1);
            app.gridvars{app.gridtab.name(1:2),'val'}=[px;py];
            app.makeTrajData();
            app.integrateTraj('go','p0')
        end
        
        %clickin anywhere else in the axis triggers click-trajectories
        function clickTraj(app,src,event)
            disp('clickTraj')
            for i=1:app.nClick
                [px,py]=ginput(1);
                coords(i,:)=[px,py];
                app.markerPquery(i).XData=px;
                app.markerPquery(i).YData=py;
                app.markerPquery(i).Visible='on';
            end
            coords
            app.makeTrajData(coords);
            app.integrateTraj('go','click');
        end
        
        
        function makeTrajData(app, coords)
            isP0=false;
            if ~exist('coords','var') %p0
                coords=app.gridtab.val';
                isP0=true;
            end
            
            newP=repmat(app.p0,size(coords,1),1);
            
            %find the coords' nearest point on the grid (for x0)
            [X,Y]=ndgrid(app.gridx, app.gridy);
            C=[X(:),Y(:)]; 
            pix = knnsearch(C,coords);
            newX0=app.XF(pix(:),:);
            
            %now overwrite relevant coords
            if app.gridtab{'x','type'}=="par"
                newP(:,app.gridtab.ix(1))=coords(:,1);
            else
                newX0(:,app.gridtab.ix(1))=coords(:,1);
            end
            if app.gridtab{'y','type'}=="par"
                newP(:,app.gridtab.ix(2))=coords(:,2);
            else
                newX0(:,app.gridtab.ix(2))=coords(:,2);
            end
            
            %store the new p0/x0: will load the appropriate one in
            %integrateTraj routine
            if isP0
                app.trajP0.p0=newP;
                app.trajP0.x0=newX0;
%                 app.gridvars{app.gridvars.type=="ic",'val'}=newX0(:);
            else
                for i=1:size(coords,1)
                    app.trajClick(i).p0=newP(i,:);
                    app.trajClick(i).x0=newX0(i,:);
                end
            end
        end
        
        
        function integrateTraj(app, action, which_traj)
            %manage tspan and x0 host side to allow two independent streams (p0 and click)
            switch which_traj
                case 'p0'
                    traj=app.trajP0;
                case 'click'
                    traj=app.trajClick;
            end
            
            tic
            append=false;
            switch action
                case 'continue'
                    if isempty(traj(1).t) %only works if prior trajectory available
                        return
                    end
                    %NOTE: non-autonomous needs time shift!
%                     tspan=app.trajTspan + traj(1).t(end);
%                     app.clo_t.settspan(tspan);
                    app.clo_t.shiftTspan();
                    for i=1:length(traj)
                        newX0(i,:)=traj(i).x(end,1:app.prob.nVar);
                    end
                    append=true;
                    
                case 'point' %no prep
                    app.clo_t.settspan(app.trajTspan);
                    for i=1:length(traj)
                        newX0(i,:)=app.gridvars.val(app.gridvars.type=="ic")';
                    end

                case 'go' %no prep
                    app.clo_t.settspan(app.trajTspan);
                    newX0 = cat(1,traj(:).x0);
                    
                case 'last'
                    if isempty(traj(1).t) %only works if prior trajectory available
                        return
                    end
                    app.clo_t.settspan(app.trajTspan);
                    for i=1:length(traj)
                        newX0(i,:)=traj(i).x(end,1:app.prob.nVar);
                    end
                    
                case 'shift'
                    if isempty(traj(1).t) %only works if prior trajectory available
                        return
                    end
%                     tspan=app.trajTspan + traj(1).t(end)
%                     app.clo_t.settspan(tspan);
                    app.clo_t.shiftTspan();
                    for i=1:length(traj)
                        newX0(i,:)=traj(i).x(end,1:app.prob.nVar);
                    end
                    
                case 'random'
                    x0lb=app.gridvars.lb(app.gridvars.type=="ic")';
                    x0ub=app.gridvars.ub(app.gridvars.type=="ic")';
                    app.clo_t.settspan(app.trajTspan);
                    for i=1:length(traj)
                        newX0(i,:)=x0lb+rand(1, length(x0lb)).*(x0ub-x0lb);
                    end
                    
%                 case 'transient'
%                     traj.x0=[traj(:).x(end,:)];
%                     app.clo_t.transient(); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            newP=cat(1,traj(:).p0);
            app.clo_t.setProblemData(newX0, newP);
            app.clo_t.trajectory(); %run the trajectory
            toc
            
            %extract results
            tt=app.clo_t.getT()/app.clo_t.tscale;
            xx=app.clo_t.getX();
            aux=app.clo_t.getAux();
            dx=app.clo_t.getDx();
            nStored=app.clo_t.getNstored();
            
            nTraj=length(traj);
            for i=1:nTraj
                traj(i).p0=newP(i,:);
                traj(i).x0=newX0(i,:);
                
                newt=tt(1:nStored(i),i);
                newx=[xx(1:nStored(i),:,i),...
                    aux(1:nStored(i),:,i),...
                    dx(1:nStored(i),:,i)];
                if append
                    newt=[traj(i).t; newt-(newt(1)-traj(i).t(end))];
                    newx=[traj(i).x; newx];
                end
                traj(i).t=newt;
                traj(i).x=newx;
            end
            
            %store the result and update plot
            switch which_traj
                case 'p0'
                    app.trajP0=traj;
                case 'click'
                    app.trajClick=traj;
            end
            
            app.updateTrajPlotData(which_traj);
            app.updateTrajPlot(which_traj);
        end
        
        function updateTrajPlotData(app, which_traj)
            switch which_traj
                case 'p0'
                    traj=app.trajP0;
                case 'click'
                    traj=app.trajClick;
            end
            if isempty(traj(1).t)
                return
            end
            nTraj=length(traj)
            for i=1:nTraj
                %extract the plotting data
                xname=app.trajXDropDown.Value;
                if xname=="t"
                    xp=traj(i).t/app.tscale;
                else
                    xix=strcmp(app.trajvars,xname);
                    xp=traj(i).x(:,xix);
                end

                yname=app.trajYDropDown.Value;
                yix=strcmp(app.trajvars,yname);
                yp=traj(i).x(:,yix);
                
                traj(i).xp=xp;
                traj(i).yp=yp;
                
                if app.trajZDropDown.Value~="none"
                    zname=app.trajZDropDown.Value;
                    zix=strcmp(app.trajvars,zname);
                    zp=traj(i).x(:,zix);
                    traj(i).zp=zp;
                end
            end
            
            %finally, store the result and update plot
            switch which_traj
                case 'p0'
                    app.trajP0=traj;
                case 'click'
                    app.trajClick=traj;
            end
        end
        
        function updateTrajPlot(app, which_traj)
            switch which_traj
                case 'p0'
                    traj=app.trajP0;
                    if ~isvalid(app.figTrajP0)
                        app.createTrajP0Fig();
                    end
                    fig=app.figTrajP0;
                    axyy=app.axyyTrajP0;
                    lyy=app.lineyyTrajP0;
                    ax3d=app.ax3DTrajP0;
                    l3d=app.line3DTrajP0;
                    
                    app.showTrajP0CheckBox.Value=1;
                    app.showTrajP0CheckBoxValueChanged()
                    
                case 'click'
                    traj=app.trajClick;
                    if ~isvalid(app.figTrajClick)
                        app.createTrajClickFig();
                    end
                    fig=app.figTrajClick;
                    axyy=app.axyyTrajClick;
                    lyy=app.lineyyTrajClick;
                    ax3d=app.ax3DTrajClick;
                    l3d=app.line3DTrajClick;
                    
                    app.showTrajClickCheckBox.Value=1;
                    app.showTrajClickCheckBoxValueChanged()
            end
            
            if isempty(traj(1).t)
                return
            end
            
            xname=app.trajXDropDown.Value;
            yname=app.trajYDropDown.Value;
            
            for i=1:length(traj)
                refreshdata(lyy(i,1))
                if app.trajZDropDown.Value~="none"
                    
                    zname=app.trajZDropDown.Value;
                    
                    refreshdata(lyy(i,2))
                    yyaxis(axyy(i),'right');
                    ylabel(axyy(i), zname);
                    
                    refreshdata(l3d(i))
                    xlabel(ax3d(i), xname);
                    ylabel(ax3d(i), yname);
                    zlabel(ax3d(i), zname);
                    
                    if app.threeDButton.Value==1
                        axyy(i).Visible='off';
                        axyy(i).YAxis(2).Visible='off'; %off until y2 is set
                        lyy(i,1).Visible='off';
                        lyy(i,2).Visible='off';
                        ax3d(i).Visible='on';
                        l3d(i).Visible='on';
                    else
                        axyy(i).Visible='on';
                        axyy(i).YAxis(2).Visible='on'; %off until y2 is set
                        lyy(i,1).Visible='on';
                        lyy(i,2).Visible='on';
                        ax3d(i).Visible='off';
                        l3d(i).Visible='off';
                    end
                else
                    axyy(i).YAxis(2).Visible='off'; %off until y2 is set
                    lyy(i,2).Visible='off';
                    ax3d(i).Visible='off';
                    l3d(i).Visible='off';
                end
                
                yyaxis(axyy(i),'left'); %bring focus back to left
                ylabel(axyy(i), yname);
            end
            xlabel(axyy(end), xname);
            if app.linkAxesButton.Value
                app.setLinkedLims();
            end
            figure(fig)
        end

    end

    
    
    % Component initialization
    methods (Access = private)
        
        function createControlFig(app)
            %TODO: normalized position units, relative to screen size?
            
            % Create figControl and hide until all components are created
            app.figControl = uifigure();
            app.figControl.Visible = 'off';
            app.figControl.Position = [50 400 300 650];
            app.figControl.Name = 'controls';
            app.figControl.CloseRequestFcn = @app.figControlCloseRequest;
            
            % Create TabGroup
            app.TabGroup = uitabgroup(app.figControl);
            app.TabGroup.Position = [5 5 290 640];
            
            % Create NumericsTab
            app.NumericsTab = uitab(app.TabGroup);
            app.NumericsTab.Title = 'Numerics';
            
            % Create SelectODEfileButton
            app.SelectODEfileButton = uibutton(app.NumericsTab, 'push');
            app.SelectODEfileButton.ButtonPushedFcn = @app.SelectODEfileButtonPushed;
            app.SelectODEfileButton.Position = [5 590 115 20];
            app.SelectODEfileButton.Text = 'Select ODE file';
            
            % Create currentFileLabel
            app.currentFileLabel = uilabel(app.NumericsTab);
            app.currentFileLabel.Position = [125 590 160 20];
            app.currentFileLabel.Text = 'currentFile';
            app.currentFileLabel.HorizontalAlignment='center';
            app.currentFileLabel.BackgroundColor='w';
            
            % Create gridDeviceDropDownLabel
            app.gridDeviceDropDownLabel = uilabel(app.NumericsTab);
            app.gridDeviceDropDownLabel.Position = [5 560 75 20];
            app.gridDeviceDropDownLabel.Text = 'Grid device';
            
            % Create gridDeviceDropDown
            app.gridDeviceDropDown = uidropdown(app.NumericsTab);
            app.gridDeviceDropDown.ValueChangedFcn = @app.gridDeviceDropDownValueChanged;
            app.gridDeviceDropDown.Position = [5 540 135 20];
            
            % Create gridPrecisionCheckBox
            app.gridPrecisionCheckBox = uicheckbox(app.NumericsTab);
            app.gridPrecisionCheckBox.ValueChangedFcn = @app.gridPrecisionCheckBoxValueChanged;
            app.gridPrecisionCheckBox.Text = 'double precision';
            app.gridPrecisionCheckBox.Position = [10 520 125 20];
            
            % Create gridStepperDropDownLabel
            app.gridStepperDropDownLabel = uilabel(app.NumericsTab);
            app.gridStepperDropDownLabel.Position = [5 500 50 20];
            app.gridStepperDropDownLabel.Text = 'stepper';
            
            % Create gridStepperDropDown
            app.gridStepperDropDown = uidropdown(app.NumericsTab);
            app.gridStepperDropDown.ValueChangedFcn = @app.gridStepperDropDownValueChanged;
            app.gridStepperDropDown.Position = [60 500 80 20];
            
            % Create gridObserverDropDownLabel
            app.gridObserverDropDownLabel = uilabel(app.NumericsTab);
            app.gridObserverDropDownLabel.Position = [5 465 55 20];
            app.gridObserverDropDownLabel.Text = 'observer';
            
            % Create gridObserverDropDown
            app.gridObserverDropDown = uidropdown(app.NumericsTab);
            app.gridObserverDropDown.ValueChangedFcn = @app.gridObserverDropDownValueChanged;
            app.gridObserverDropDown.Position = [60 465 80 20];
            
            % Create gridCLIsBuiltButton
            app.gridCLIsBuiltButton = uibutton(app.NumericsTab, 'state');
            app.gridCLIsBuiltButton.ValueChangedFcn = @app.gridCLIsBuiltButtonValueChanged;
            app.gridCLIsBuiltButton.Text = 'Build CL';
            app.gridCLIsBuiltButton.Value=0;
            app.gridCLIsBuiltButton.BackgroundColor=[1,.7,.7];
            app.gridCLIsBuiltButton.Position = [30 430 85 25];
            
            
            % Create trajDeviceDropDownLabel
            app.trajDeviceDropDownLabel = uilabel(app.NumericsTab);
            app.trajDeviceDropDownLabel.Position = [150 560 140 20];
            app.trajDeviceDropDownLabel.Text = 'Trajectory device';
            
            % Create trajDeviceDropDown
            app.trajDeviceDropDown = uidropdown(app.NumericsTab);
            app.trajDeviceDropDown.ValueChangedFcn = @app.trajDeviceDropDownValueChanged;
            app.trajDeviceDropDown.Position = [150 540 135 20];
            
            % Create trajPrecisionCheckBox
            app.trajPrecisionCheckBox = uicheckbox(app.NumericsTab);
            app.trajPrecisionCheckBox.ValueChangedFcn = @app.trajPrecisionCheckBoxValueChanged;
            app.trajPrecisionCheckBox.Text = 'double precision';
            app.trajPrecisionCheckBox.Position = [155 520 125 20];
            
            % Create trajStepperDropDownLabel
            app.trajStepperDropDownLabel = uilabel(app.NumericsTab);
            app.trajStepperDropDownLabel.HorizontalAlignment = 'center';
            app.trajStepperDropDownLabel.Position = [150 500 50 20];
            app.trajStepperDropDownLabel.Text = 'stepper';
            
            % Create trajStepperDropDown
            app.trajStepperDropDown = uidropdown(app.NumericsTab);
            app.trajStepperDropDown.ValueChangedFcn = @app.trajStepperDropDownValueChanged;
            app.trajStepperDropDown.Position = [205 500 80 20];
            
            % Create trajCLIsBuiltButton
            app.trajCLIsBuiltButton = uibutton(app.NumericsTab, 'state');
            app.trajCLIsBuiltButton.ValueChangedFcn = @app.trajCLIsBuiltButtonValueChanged;
            app.trajCLIsBuiltButton.Text = 'Build CL';
            app.trajCLIsBuiltButton.Value=0;
            app.trajCLIsBuiltButton.BackgroundColor=[1,.7,.7];
            app.trajCLIsBuiltButton.Position = [175 430 85 25];
            
            
            
            % Create SolverParLabel
            app.SolverParLabel = uilabel(app.NumericsTab);
            app.SolverParLabel.Position = [5 400 150 20];
            app.SolverParLabel.Text = 'Solver parameters';
            
            % Create SolverParTable
            app.SolverParTable = uitable(app.NumericsTab);
            app.SolverParTable.ColumnName = {'grid','traj'};
            app.SolverParTable.ColumnWidth = '1x';
            app.SolverParTable.ColumnEditable = [true, true];
            app.SolverParTable.CellEditCallback = @app.SolverParTableCellEdit;
            app.SolverParTable.Position = [5 215 280 185];
            
            
            % Create gridObserverParLabel
            app.gridObserverParLabel = uilabel(app.NumericsTab);
            app.gridObserverParLabel.Position = [5 190 150 20];
            app.gridObserverParLabel.Text = 'Observer parameters';
            
            % Create gridObserverParTable
            app.gridObserverParTable = uitable(app.NumericsTab);
            app.gridObserverParTable.ColumnName = {'value'};
            app.gridObserverParTable.ColumnWidth = '1x';
            app.gridObserverParTable.ColumnEditable = [true];
            app.gridObserverParTable.CellEditCallback = @app.gridObserverParTableCellEdit;
            app.gridObserverParTable.Position = [5 5 280 185];
            
            
            
            % Create GridTab
            app.GridTab = uitab(app.TabGroup);
            app.GridTab.Title = 'Grid';
            
            % Create gridTable
            app.gridTable = uitable(app.GridTab);
            app.gridTable.RowName={'x';'y'};
            app.gridTable.ColumnName={'name','val','lb','ub','N'}; %,'del'
            app.gridTable.ColumnWidth = '1x';
            app.gridTable.ColumnEditable = [true true true true true];
            app.gridTable.CellEditCallback = @app.gridTableCellEdit;
            app.gridTable.Position = [5 535 280 75];
            
            % Create gridZDropDownLabel
            app.gridZDropDownLabel = uilabel(app.GridTab);
            app.gridZDropDownLabel.Position = [5 510 10 20];
            app.gridZDropDownLabel.Text = 'z';
            app.gridZDropDownLabel.HorizontalAlignment = 'right';
            
            % Create gridZDropDown
            app.gridZDropDown = uidropdown(app.GridTab);
            app.gridZDropDown.ValueChangedFcn = @app.gridZDropDownValueChanged;
            app.gridZDropDown.Position = [20 510 65 20];
            
            % Create zValueEditField
            app.zValueEditField = uieditfield(app.GridTab, 'numeric');
            app.zValueEditField.ValueChangedFcn = @app.zValueEditFieldValueChanged;
            app.zValueEditField.Position = [90 510 60 20];
            app.zValueEditField.Value=0;
            app.zValueEditField.ValueDisplayFormat='%.2g';
            
            % Create dzEditFieldLabel
            app.dzEditFieldLabel = uilabel(app.GridTab);
            app.dzEditFieldLabel.Position = [155 510 15 20];
            app.dzEditFieldLabel.Text = 'dz';
            app.dzEditFieldLabel.HorizontalAlignment = 'left';
            
            % Create dzEditField
            app.dzEditField = uieditfield(app.GridTab, 'numeric');
            app.dzEditField.ValueChangedFcn = @app.dzEditFieldValueChanged;
            app.dzEditField.Position = [170 510 60 20];
            app.dzEditField.Value=0;
            app.dzEditField.ValueDisplayFormat='%.2g';
            
            % Create featureDropDownLabel
            app.featureDropDownLabel = uilabel(app.GridTab);
            app.featureDropDownLabel.Position = [200 490 85 20];
            app.featureDropDownLabel.Text = 'feature';
            app.featureDropDownLabel.HorizontalAlignment='right';
            
            % Create featureDropDown
            app.featureDropDown = uidropdown(app.GridTab);
            app.featureDropDown.ValueChangedFcn = @app.featureDropDownValueChanged;
            app.featureDropDown.Position = [160 474 125 20];
            
            % Create fscaleEditFieldLabel
            app.fscaleEditFieldLabel = uilabel(app.GridTab);
            app.fscaleEditFieldLabel.Position = [195 453 40 20];
            app.fscaleEditFieldLabel.Text = 'fscale';
            
            % Create fscaleEditField
            app.fscaleEditField = uieditfield(app.GridTab, 'numeric');
            app.fscaleEditField.ValueChangedFcn = @app.fscaleEditFieldValueChanged;
            app.fscaleEditField.Position = [235 453 50 20];
            app.fscaleEditField.Value=1;
            
            % Create tscaleEditFieldLabel
            app.tscaleEditFieldLabel = uilabel(app.GridTab);
            app.tscaleEditFieldLabel.Position = [195 431 40 20];
            app.tscaleEditFieldLabel.Text = 'tscale';
            
            % Create tscaleEditField
            app.tscaleEditField = uieditfield(app.GridTab, 'numeric');
            app.tscaleEditField.ValueChangedFcn = @app.tscaleEditFieldValueChanged;
            app.tscaleEditField.Position = [235 431 50 20];
            app.tscaleEditField.Value=1;
            
            % Create tspanTable
            app.tspanTable = uitable(app.GridTab);
            app.tspanTable.ColumnName = {'t0','tf'};
            app.tspanTable.RowName = {'grid','traj'};
            app.tspanTable.ColumnWidth = '1x';
            app.tspanTable.ColumnEditable = [true true];
            app.tspanTable.CellEditCallback = @app.tspanTableCellEdit;
            app.tspanTable.Position = [5 430 150 75];
            
            
            % Create parLabel
            app.parLabel = uilabel(app.GridTab);
            app.parLabel.Position = [5 400 69 23];
            app.parLabel.Text = 'Parameters';
            
            % Create parSaveButton
            app.parSaveButton = uibutton(app.GridTab, 'push');
            app.parSaveButton.ButtonPushedFcn = @app.parSaveButtonPushed;
            app.parSaveButton.Position = [180 403 48 20];
            app.parSaveButton.Text = 'save';
            
            % Create parDefaultButton
            app.parDefaultButton = uibutton(app.GridTab, 'push');
            app.parDefaultButton.ButtonPushedFcn = @app.parDefaultButtonPushed;
            app.parDefaultButton.Position = [233 403 52 20];
            app.parDefaultButton.Text = 'default';
            
            % Create parTable
            app.parTable = uitable(app.GridTab);
            app.parTable.RowName={};
            app.parTable.ColumnName = {'name', 'val', 'lb', 'ub'};
            app.parTable.ColumnWidth = '1x';
            app.parTable.ColumnSortable = [true false false false];
            app.parTable.ColumnEditable = [true true true true];
            app.parTable.CellEditCallback = @app.parTableCellEdit;
            app.parTable.Position = [5 225 280 175];
            
            
            % Create icLabel
            app.icLabel = uilabel(app.GridTab);
            app.icLabel.Position = [5 195 91 23];
            app.icLabel.Text = 'Initial conditions';
            
            % Create icDefaultButton
            app.icDefaultButton = uibutton(app.GridTab, 'push');
            app.icDefaultButton.ButtonPushedFcn = @app.icDefaultButtonPushed;
            app.icDefaultButton.Position = [233 198 52 20];
            app.icDefaultButton.Text = 'default';
            
            % Create icTable
            app.icTable = uitable(app.GridTab);
            app.icTable.RowName={};
            app.icTable.ColumnName = {'name', 'val', 'lb', 'ub'};
            app.icTable.ColumnWidth = '1x';
            app.icTable.ColumnSortable = [true false false false];
            app.icTable.ColumnEditable = [true true true true];
            app.icTable.CellEditCallback = @app.icTableCellEdit;
            app.icTable.Position = [5 75 280 120];
            
            
            %             % Create icRandomButton
            %             app.icRandomButton = uibutton(app.GridTab, 'push');
            %             app.icRandomButton.ButtonPushedFcn = @app.icRandomButtonPushed;
            %             app.icRandomButton.Position = [229 172 52 22];
            %             app.icRandomButton.Text = 'random';
            
            
            %             app.TrajTab = uitab(app.TabGroup);
            %             app.TrajTab.Title = 'Trajectories';
            
            % Create trajectoriesLabel
            app.trajectoriesLabel = uilabel(app.GridTab);
            app.trajectoriesLabel.Position = [5 50 68 23];
            app.trajectoriesLabel.Text = 'Trajectories';
            
            % Create showTrajP0CheckBox
            app.showTrajP0CheckBox = uicheckbox(app.GridTab);
            app.showTrajP0CheckBox.ValueChangedFcn = @app.showTrajP0CheckBoxValueChanged;
            app.showTrajP0CheckBox.Text = 'p0';
            app.showTrajP0CheckBox.Position = [5 30 55 22];
            
            % Create showTrajClickCheckBox
            app.showTrajClickCheckBox = uicheckbox(app.GridTab);
            app.showTrajClickCheckBox.ValueChangedFcn = @app.showTrajClickCheckBoxValueChanged;
            app.showTrajClickCheckBox.Text = 'click';
            app.showTrajClickCheckBox.Position = [60 30 55 22];
            
            % Create clicksEditFieldLabel
            app.clicksEditFieldLabel = uilabel(app.GridTab);
            app.clicksEditFieldLabel.HorizontalAlignment = 'center';
            app.clicksEditFieldLabel.Position = [115 30 45 22];
            app.clicksEditFieldLabel.Text = '# clicks';
            
            % Create clicksEditField
            app.clicksEditField = uieditfield(app.GridTab, 'numeric');
            app.clicksEditField.ValueChangedFcn = @app.clicksEditFieldValueChanged;
            app.clicksEditField.Position = [160 30 30 22];
            app.clicksEditField.Value=app.nClick;
            
            % Create linkAxesButton
            app.linkAxesButton = uibutton(app.GridTab, 'state');
            app.linkAxesButton.ValueChangedFcn = @app.linkAxesButtonValueChanged;
            app.linkAxesButton.Text = 'link axes';
            app.linkAxesButton.Position = [195 30 55 22];
            
            % Create threeDButton
            app.threeDButton = uibutton(app.GridTab, 'state');
            app.threeDButton.ValueChangedFcn = @app.threeDButtonValueChanged;
            app.threeDButton.Text = '3D';
            app.threeDButton.Position = [255 30 30 22];
            
            % Create trajXDropDownLabel
            app.trajXDropDownLabel = uilabel(app.GridTab);
            app.trajXDropDownLabel.HorizontalAlignment = 'center';
            app.trajXDropDownLabel.Position = [5 5 15 20];
            app.trajXDropDownLabel.Text = 'x';
            
            % Create trajXDropDown
            app.trajXDropDown = uidropdown(app.GridTab);
            app.trajXDropDown.ValueChangedFcn = @app.trajXDropDownValueChanged;
            app.trajXDropDown.Position = [20 5 75 20];
            
            % Create trajYDropDownLabel
            app.trajYDropDownLabel = uilabel(app.GridTab);
            app.trajYDropDownLabel.HorizontalAlignment = 'left';
            app.trajYDropDownLabel.Position = [100 5 10 20];
            app.trajYDropDownLabel.Text = 'y1';
            
            % Create trajYDropDown
            app.trajYDropDown = uidropdown(app.GridTab);
            app.trajYDropDown.ValueChangedFcn = @app.trajYDropDownValueChanged;
            app.trajYDropDown.Position = [115 5 75 20];
            
            % Create trajZDropDownLabel
            app.trajZDropDownLabel = uilabel(app.GridTab);
            app.trajZDropDownLabel.HorizontalAlignment = 'left';
            app.trajZDropDownLabel.Position = [195 5 15 20];
            app.trajZDropDownLabel.Text = 'y2';
            
            % Create trajZDropDown
            app.trajZDropDown = uidropdown(app.GridTab);
            app.trajZDropDown.ValueChangedFcn = @app.trajZDropDownValueChanged;
            app.trajZDropDown.Position = [210 5 75 20];
            
            %             app.HelpTab = uitab(app.TabGroup);
            %             app.HelpTab.Title = 'Usage';
            
            app.figControl.KeyPressFcn = @app.commonKeyPress;

            app.figControl.Visible = 'on';
            
        end
        
        function createGridFig(app)
            % Create figGrid and hide until all components are created
            app.figGrid = figure(1983);
            app.figGrid.Visible = 'off';
            app.figGrid.Position = [375 400 600 600];
            app.figGrid.Name = 'grid';
            app.figGrid.NumberTitle='off';
            
            % grid axis
            t=tiledlayout(app.figGrid,1,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axGrid = nexttile(t);
            
            app.imGrid=imagesc(app.axGrid,0, 0, 0);
            app.imGrid.HitTest='off';  %not clickable
            
            hold on
            [~,app.gridContours] = contour(app.axGrid, [0,1], [0,1], [0,0;0,1],ShowText="on");
            app.gridContours.EdgeColor=0.5*[1,1,1];
            app.gridContours.Visible="off";
            hold off

            app.axGrid.YDir='normal';
            axis(app.axGrid,'tight');
            app.gridCBar=colorbar('northoutside');
            colormap(app.axGrid, app.gridColormap)
%             colormap(app.axGrid,[0.85*[1,1,1];turbo])
            
            xlabel(app.axGrid,'x')
            ylabel(app.axGrid,'y')
            
            app.axGrid.ButtonDownFcn = @app.clickTraj; %axis is clickable
            app.figGrid.KeyPressFcn = @app.commonKeyPress;
            app.figGrid.Visible = 'on';
            grid(app.axGrid,"on")
        end
        
        function createTrajP0Fig(app)
            % Create figTraj and hide until all components are created
            app.figTrajP0 = figure(1984);
            app.figTrajP0.Visible = 'off';
            app.figTrajP0.Position = [375 100 600 200];
            app.figTrajP0.Name = 'p0';
            app.figTrajP0.NumberTitle='off';
            

            t=tiledlayout(app.figTrajP0,1,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axyyTrajP0 = nexttile(t);
            app.ax3DTrajP0 = copyobj(app.axyyTrajP0,app.figTrajP0);
            
            yyaxis(app.axyyTrajP0,'left');
            app.lineyyTrajP0(1)=plot(app.axyyTrajP0,nan,nan);
            yyaxis(app.axyyTrajP0,'right');
            app.lineyyTrajP0(2)=plot(app.axyyTrajP0,nan,nan);
            
            app.lineyyTrajP0(1).XDataSource='app.trajP0.xp';
            app.lineyyTrajP0(1).YDataSource='app.trajP0.yp';
            app.lineyyTrajP0(2).XDataSource='app.trajP0.xp';
            app.lineyyTrajP0(2).YDataSource='app.trajP0.zp';
            
            xlabel(app.axyyTrajP0,'t')
            axis(app.axyyTrajP0,'tight')
            app.axyyTrajP0.YAxis(2).Visible='off'; %turn it on if y2 is set
            
            app.line3DTrajP0=plot3(app.ax3DTrajP0,nan,nan,nan);
            app.line3DTrajP0.XDataSource='app.trajP0.xp';
            app.line3DTrajP0.YDataSource='app.trajP0.yp';
            app.line3DTrajP0.ZDataSource='app.trajP0.zp';
            axis(app.ax3DTrajP0,'vis3d')
            
            app.ax3DTrajP0.Visible='off'; %start in yyaxis mode
            app.figTrajP0.KeyPressFcn = @app.commonKeyPress;
            app.figTrajP0.Visible = 'on';

            app.markerP0=line(app.axGrid,nan,nan,'color','k','marker','o','linestyle','none');
            app.markerP0.ButtonDownFcn=@app.clickP0;  %p0 marker is clickable
        end
        
        function createTrajClickFig(app)
            % Create figTraj and hide until all components are created
            app.figTrajClick = figure(1985);
            app.figTrajClick.Visible = 'off';
            app.figTrajClick.Position = [1000 100 900 900];
            app.figTrajClick.Name = 'click trajectories';
            app.figTrajClick.NumberTitle='off';
            
            app.tilesTrajClick=tiledlayout(app.figTrajClick,app.nClick,1);
            app.tilesTrajClick.TileSpacing = 'compact';
            app.tilesTrajClick.Padding = 'compact';
            
%             app.axyyTrajClick(app.nClick,1)=matlab.graphics.axis.Axes;
%             app.lineyyTrajClick(app.nClick,2)=matlab.graphics.chart.primitive.Line;
%             app.ax3DTrajClick(app.nClick,1)=matlab.graphics.axis.Axes;
%             app.line3DTrajClick(app.nClick,1)=matlab.graphics.chart.primitive.Line;
            
            for i=1:app.nClick
                app.axyyTrajClick(i,1) = nexttile(app.tilesTrajClick);
                app.ax3DTrajClick(i,1) = copyobj(app.axyyTrajClick(i),app.figTrajClick);
                app.ax3DTrajClick(i,1).Position=app.axyyTrajClick(i,1).Position;
                
                yyaxis(app.axyyTrajClick(i),'left');
                app.lineyyTrajClick(i,1)=plot(app.axyyTrajClick(i),nan,nan);
                yyaxis(app.axyyTrajClick(i),'right');
                app.lineyyTrajClick(i,2)=plot(app.axyyTrajClick(i),nan,nan);
                
                app.lineyyTrajClick(i,1).XDataSource=['app.trajClick(',num2str(i),').xp'];
                app.lineyyTrajClick(i,1).YDataSource=['app.trajClick(',num2str(i),').yp'];
                app.lineyyTrajClick(i,2).XDataSource=['app.trajClick(',num2str(i),').xp'];
                app.lineyyTrajClick(i,2).YDataSource=['app.trajClick(',num2str(i),').zp'];

                app.axyyTrajClick(i).XTickLabel=[];
                axis(app.axyyTrajClick(i),'tight')
%                 app.axyyTrajClick(i).YAxis(2).Visible='off'; %off until y2 is set
                yyaxis(app.axyyTrajClick(i),'left'); %focus on left
                
                app.line3DTrajClick(i)=plot3(app.ax3DTrajClick(i),nan,nan,nan);
                app.line3DTrajClick(i).XDataSource=['app.trajClick(',num2str(i),').xp'];
                app.line3DTrajClick(i).YDataSource=['app.trajClick(',num2str(i),').yp'];
                app.line3DTrajClick(i).ZDataSource=['app.trajClick(',num2str(i),').zp'];
                axis(app.ax3DTrajClick(i),'tight')
                app.ax3DTrajClick(i).Visible='off'; %start in yyaxis mode
                axis(app.ax3DTrajClick(i),'vis3d');
            end 
            app.axyyTrajClick(end).XTickLabelMode='auto';
            
            app.figTrajClick.KeyPressFcn = @app.commonKeyPress;
            app.figTrajClick.Visible = 'on';

            markers='sdp^v';
            for i=1:app.nClick
                app.markerPquery(i)=line(app.axGrid,nan,nan,'color','k','marker',markers(i),'linestyle','none');
                app.markerPquery(i).HitTest='off'; %not clickable
            end
        end
    end
    
end