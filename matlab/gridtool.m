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
        z %struct with z info: row from gridvars
        dz %increment to change z when using +/-
        fscale
        XF %place to store the final state of most recent grid simulation
        
        nClick=3
        
        %mechanism to cache F(:,:,z) values - for rapidly scanning z post compute
        
    end
    
    properties (Dependent = true)
        gridx %grid xcoords
        gridy %grid ycoords
        p0 %current parameter vector value (extract from gridvars)
        x0 %current initial condition value
        P0 %base parameter matrix with p0 replicated nPts times [nPts, nPar]
        X0 %base ic matrix with x0 replicated nPts times
        nPts
    end
    
    %data needed to display graphs: set-observable, so that whenever these
    %change, triggers update to relevant figures
    properties (SetObservable = true)
        
        gridvars %table of possible grid vars (par+var) with values
        grid=table('Size',[2,6], 'VariableTypes',...
            {'cellstr','double','double','double','categorical','double'},...
            'VariableNames',{'name','val','lb','ub','type','ix'},...
            'RowNames',{'x','y'}) %table specifying grid x, y info
        feature %struct specifying feature info for grid display
        
%         gridSol=struct('x',[],'y',[],'F',[]);
        trajSol=struct('t',[],'x',[],'xname',{},'tlim',[],'xlo',[],'xhi',[]); %struct for trajectory data
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
        figGrid
        axGrid                  matlab.graphics.axis.Axes
        imGrid                  matlab.graphics.primitive.Image
        gridCBar
        markerP0 %marker (line object) specifying (x,y) for p0
        markerPquery %markers (line objects) specifying (x,y) for n query trajectories
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %figure with p0 trajectory
        figTrajP0
        axyyTrajP0                matlab.graphics.axis.Axes
        lineyyTrajP0              matlab.graphics.chart.primitive.Line 
        ax3DTrajP0                matlab.graphics.axis.Axes
        line3DTrajP0              matlab.graphics.chart.primitive.Line 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %click trajectory figure with n trajectories
        figTrajClick
        tilesTrajClick
        axyyTrajClick             matlab.graphics.axis.Axes
        lineyyTrajClick           matlab.graphics.chart.primitive.Line 
        ax3DTrajClick             matlab.graphics.axis.Axes
        line3DTrajClick           matlab.graphics.chart.primitive.Line 
        
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
        
        SolverParLabel      matlab.ui.control.Label
        SolverParTable      matlab.ui.control.Table
        gridObserverParTable    matlab.ui.control.Table
        gridObserverParLabel    matlab.ui.control.Label
        
        GridTab                 matlab.ui.container.Tab
        
        gridTable               matlab.ui.control.Table
        gridZDropDownLabel      matlab.ui.control.Label
        gridZDropDown           matlab.ui.control.DropDown
        dzEditFieldLabel        matlab.ui.control.Label
        dzEditField             matlab.ui.control.NumericEditField
        tspanTable              matlab.ui.control.Table
        featureDropDownLabel    matlab.ui.control.Label
        featureDropDown         matlab.ui.control.DropDown
        fscaleEditFieldLabel    matlab.ui.control.Label
        fscaleEditField         matlab.ui.control.NumericEditField
        parLabel                matlab.ui.control.Label
        parTable                matlab.ui.control.Table
        icLabel                 matlab.ui.control.Label
        icTable                 matlab.ui.control.Table
        parDefaultButton        matlab.ui.control.Button
        icDefaultButton         matlab.ui.control.Button
        
%         TrajTab                 matlab.ui.container.Tab
        
        trajectoriesLabel       matlab.ui.control.Label
        clicksEditFieldLabel    matlab.ui.control.Label
        clicksEditField         matlab.ui.control.NumericEditField
        trajXDropDownLabel          matlab.ui.control.Label
        trajXDropDown               matlab.ui.control.DropDown
        trajYDropDownLabel          matlab.ui.control.Label
        trajYDropDown               matlab.ui.control.DropDown
        trajZDropDownLabel          matlab.ui.control.Label
        trajZDropDown               matlab.ui.control.DropDown
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
        
        % Value changed function: gridCLIsBuiltButton (doesn't execute for
        % programmatic changes)
        function gridCLIsBuiltButtonValueChanged(app, src, event)
            isBuilt = app.gridCLIsBuiltButton.Value;
            if ~isBuilt
                app.buildCL('grid');
                app.clo_g.initialize();
            end
        end
        
        % Value changed function: trajCLIsBuiltButton (doesn't execute for
        % programmatic changes)
        function trajCLIsBuiltButtonValueChanged(app, src, event)
            isBuilt = app.trajCLIsBuiltButton.Value;
            if ~isBuilt
                app.buildCL('traj');
                app.clo_t.initialize();
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
        
        
        % Cell edit callback: gridTable
        function gridTableCellEdit(app, src, event)
            % gridTable val/lb/ub must be synced to parTable and icTable
            indices = event.Indices;
            newData = event.NewData;
            vix=indices(1);
            prop=src.ColumnName{indices(2)};
            switch prop
                case 'name'
                    app.grid(vix,:)=app.gridvars(newData,:);
                    
                case {'val', 'lb', 'ub'}
                    app.gridvars{app.grid.name(vix),prop}=newData; %sync gridvars table
                    
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
        end
        
        % Value changed function: dzEditField
        function dzEditFieldValueChanged(app, src, event)
            %TODO: value checking
            newDZ = app.dzEditField.Value;
            if isnumeric(newDZ) && isreal(newDZ)
                app.dz=newDZ;
            end
        end
        
        
        % Cell edit callback: tspanTable
        function tspanTableCellEdit(app, src, event)
            indices = event.Indices; %r=1:grid, r=2:traj, c=1:t0, c2:tf, c3:trans
            newData = event.NewData;
            if indices(1)==1
                tspan=app.clo_g.tspan;
                tspan(indices(2))=newData;
                app.clo_g.settspan(tspan);
            elseif indices(1)==2
                tspan=app.clo_t.tspan;
                tspan(indices(2))=newData;
                app.clo_t.settspan(tspan);
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
        
        % Button pushed function: parDefaultButton
        function parDefaultButtonPushed(app, src, event)
            %get the original default values for par.val/lb/ub from app.prob
            newGridvars=table;
            newGridvars.val=app.prob.p0(:);
            newGridvars.lb=[app.prob.par.lb]';
            newGridvars.ub=[app.prob.par.ub]';
            app.gridvars(app.gridvars.type=="par",{'val','lb','ub'})=newGridvars;
            app.makeGridData(); %in case grid wasn't changed, but other pars were
        end
        
        % Button pushed function: icDefaultButton
        function icDefaultButtonPushed(app, src, event)
            %get the original default values for var.val/lb/ub from app.prob
            newGridvars=table;
            newGridvars.val=app.prob.x0(:);
            newGridvars.lb=[app.prob.var.lb]';
            newGridvars.ub=[app.prob.var.ub]';
            app.gridvars(app.gridvars.type=="ic",{'val','lb','ub'})=newGridvars;
        end
        
        % Value changed function: trajXDropDown
        function trajXDropDownValueChanged(app, src, event)
            newX = app.trajXDropDown.Value;
            
        end
        
        % Value changed function: trajYDropDown
        function trajYDropDownValueChanged(app, src, event)
            newY = app.trajYDropDown.Value;
            
        end
        
        % Value changed function: trajZDropDown
        function trajZDropDownValueChanged(app, src, event)
            newZ = app.trajZDropDown.Value;
            
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
            gridx=linspace(app.grid.lb(1),app.grid.ub(1),app.nGrid(1));
        end
        
        function gridy=get.gridy(app)
            gridy=linspace(app.grid.lb(2),app.grid.ub(2),app.nGrid(2));
        end
    end
    
    %listeners
    methods
        function attachListeners(app)
            app.listenerGrid=addlistener(app,'grid','PostSet',@gridtool.gridUpdate);
            app.listenerGridVars=addlistener(app,'gridvars','PostSet',@gridtool.gridvarUpdate);
        end
    end
    
    methods (Static = true)
        function gridUpdate(~,eventData)
            app = eventData.AffectedObject;
            newGridDisplayData=app.grid(:,{'name','val','lb','ub'});
            newGridDisplayData.N=app.nGrid(:);
            app.gridTable.Data=table2cell(newGridDisplayData);
            app.updateGridPlot('nosol');
            app.makeGridData(); %prep for next integration
        end
        
        function gridvarUpdate(~,eventData)
            app = eventData.AffectedObject;
            %update grid
            if isempty(app.grid.name{1})
                vars=app.gridvars.name(1:2);
            else
                vars=app.grid.name;
            end
            newGrid=app.gridvars(vars,:);
            newGrid.Row={'x';'y'};
            if ~isequal(newGrid, app.grid)
                app.grid(:,:)=newGrid; %grid update triggers makeGridData
            end
            newZ=app.gridvars(app.z.name,:);
            if ~isequal(newZ, app.z)
                app.z=newZ; %grid update triggers makeGridData
            end
            %update parameter and ic table
            app.parTable.Data=app.gridvars(app.gridvars.type=="par",{'name','val','lb','ub'});
            app.icTable.Data=app.gridvars(app.gridvars.type=="ic",{'name','val','lb','ub'});
            
        end
    end

    % internal functions (clODE control & ui interaction)
    methods (Access = private)
        
        function gridKeyPress(app,src,event)
            %             app.clo_g.initialize();
            disp(event.Key)
            switch(event.Key)
                case 'c' %'continue' - features without initialization
                    app.integrateGrid('continue')
                    
                case 'g' %'go' - features with initialization
                    app.integrateGrid('go')
                    
                case 'r' %'random'
                    app.integrateGrid('random')
                    
                case 't' %'transient'
                    app.integrateGrid('transient')
                    
                case 'add' % increment z by +dz, transient
                    newZ=min(app.z.val+app.dz, app.z.ub);
                    if newZ~=app.z.val
                        app.z.val=newZ;
                        app.gridvars{app.z.name,'val'}=app.z.val;
                        app.makeGridData(); %prep for next integration
                        app.integrateGrid('transient')
                    end
                    
                case 'subtract' % increment z by -dz, transient
                    newZ=max(app.z.val-app.dz, app.z.lb);
                    if newZ~=app.z.val
                        app.z.val=newZ;
                        app.gridvars{app.z.name,'val'}=app.z.val;
                        app.makeGridData(); %prep for next integration
                        app.integrateGrid('transient')
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
            
            app.clo_g.setProblemData(newX0, newP);
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
                case {'continue','go'}
                    app.updateGridPlot('feature')
                case {'random','transient'}
                    app.updateGridPlot('transient')
            end
        end
        
        function updateGridPlot(app,type)
            if ~isvalid(app.figGrid)
                app.createGridFig();
                app.imGrid.XData=app.gridx;
                app.imGrid.YData=app.gridy;
                app.markerP0.XData=app.grid{'x','val'};
                app.markerP0.YData=app.grid{'y','val'};
                xlabel(app.axGrid,app.grid{'x','name'});
                ylabel(app.axGrid,app.grid{'y','name'});
                title(app.gridCBar,app.featureDropDown.Value,'Interpreter','none')
            end
            switch type
                case 'transient'
                    %handle XF plotting selection
                    ix=app.clo_g.op.fVarIx;
                    C=reshape(app.XF(:,ix),app.nGrid)';
                    app.imGrid.CData=C;
                    title(app.gridCBar,app.prob.varNames(ix),'Interpreter','none')
                    
                case 'feature'
                    %handle F plotting seletion
                    ix=strcmp(app.featureDropDown.Value,app.clo_g.fNames);
                    C=reshape(app.clo_g.F(:,ix)*app.fscale,app.nGrid)';
                    C(C<-1e10|C>1e10)=nan; %hack to prevent display of bad values
                    app.imGrid.CData=C;
                    title(app.gridCBar,app.featureDropDown.Value,'Interpreter','none')
                    
                case 'nosol'  %setup plotting elements when grid is changed
                    app.imGrid.XData=app.gridx;
                    app.imGrid.YData=app.gridy;
                    app.imGrid.CData=zeros(app.nGrid(1),app.nGrid(2));
                    app.markerP0.XData=app.grid{'x','val'};
                    app.markerP0.YData=app.grid{'y','val'};
                    xlabel(app.axGrid,app.grid{'x','name'});
                    ylabel(app.axGrid,app.grid{'y','name'});
                    title(app.gridCBar,app.featureDropDown.Value,'Interpreter','none')
            end
            
            if app.imGrid.Visible=="off"
                app.imGrid.Visible='on';
            end
        end
        
        
        function trajKeyPress(app,src,event)
            disp(event.Key)
            switch(event.Key)
                case 'c' %'continue' - append+shift
                    app.integrateTraj('continue')
                    
                case 'g' %'go' - start here
                    app.integrateTraj('go')
                    
                case 'p' %'p0' - 
                    app.makeTrajData()
                    app.integrateTraj('go')
                    
                case 'r' %'random'
                    app.integrateTraj('random')
                    
                case 's' %'shift'
                    app.integrateTraj('shift')
                    
                case 't' %'transient'
                    app.integrateTraj('transient')
            end
        end
        
        %clicking on the p0 marker lets the user select new p0 coords
        function clickP0(app,src,event)
            disp('clickP0')
            [px,py]=ginput(1);
            app.gridvars{app.grid.name(1),'val'}=px;
            app.gridvars{app.grid.name(2),'val'}=py;
%             app.makeTrajData(); %called by grid listener
        end
        
        %clickin anywhere else in the axis triggers click-trajectories
        function clickTraj(app,src,event)
            disp('clickTraj')
            for i=1:length(app.nClick)
                [px,py]=ginput(1);
                coords(i,:)=[px,py];
            end
            app.makeTrajData(coords);
        end
        
        
        function makeTrajData(app, coords)
            if ~exist('coords','var') %clickCoords
                coords=app.grid.val';
            end
            
            newP=repmat(app.p0,size(coords,1),1);
            
            %find the coords' nearest point on the grid (for x0)
            [X,Y]=ndgrid(app.gridx, app.gridy);
            C=[X(:),Y(:)]; 
            pix = knnsearch(C,coords);
            for i=1:size(coords,1)
                newX0(i,:)=app.XF(pix(i),:);
            end
            
            %now overwrite relevant coords
            if app.grid{'x','type'}=="par"
                newP(:,app.grid.ix(1))=coords(:,1);
            else
                newX0(:,app.grid.ix(2))=coords(:,1);
            end
            if app.grid{'y','type'}=="par"
                newP(:,app.grid.ix(1))=coords(:,2);
            else
                newX0(:,app.grid.ix(2))=coords(:,2);
            end
            
            app.clo_t.setProblemData(newX0, newP);
        end
        
        
        function integrateTraj(app)
            tic
            append=false;
            switch action
                case 'continue'
                    app.clo_t.shiftX0();
                    append=true;
                    
                case 'go' %no prep
                    
                case 'shift'
                    app.clo_t.shiftTspan();
                    app.clo_t.shiftX0();
                    
                case 'random'
                    x0lb=app.gridvars.lb(app.gridvars.type=="ic")';
                    x0ub=app.gridvars.ub(app.gridvars.type=="ic")';
                    newX0=x0lb+rand(size(app.clo_t.P,1), length(x0lb)).*(x0ub-x0lb);
                    app.clo_t.setX0(newX0);
                    
                case 'transient'
                    app.clo_t.shiftX0();
            end
            app.clo_t.trajectory(); %run the trajectory
            app.clo_t.getXf();
            toc
            
            %extract results
            xx=clo.getX();
            aux=clo.getAux();
            tt=clo.getT()*clo.tscale;
            nStored=clo.getNstored();
            
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
            
        end
        
        function updateTrajP0Plot(app)
        end
        
        function updateTrajClickPlot(app)
        end
        
        function processNewODEfile(app, odefile)
            
            app.odefile=odefile;
            [~,app.prob]=ode2cl(odefile);
            
            if app.prob.nPar+app.prob.nVar<2
                error('Not enough parameters and variables for a 2D grid!')
            end
            
            app.currentFileLabel.Text=app.prob.name+".ode";
            
            %If first time, set up clODE objects using their default parameters
            if isempty(app.clo_g)
                app.clo_g=clODEfeatures(app.prob, [], app.grid_device);
                app.clo_t=clODEtrajectory(app.prob, [], app.traj_device);
                
                app.gridObserverDropDown.Items=app.clo_g.observerNames;
                app.gridObserverDropDown.Value=app.clo_g.observer;
                app.gridStepperDropDown.Items=app.clo_g.getAvailableSteppers();
                app.gridStepperDropDown.Value=app.clo_g.stepper;
                app.trajStepperDropDown.Items=app.clo_t.getAvailableSteppers();
                app.trajStepperDropDown.Value=app.clo_t.stepper;
                
                defaultsp=clODE.defaultSolverParams(); %app.clo_g.sp
                spvalues=cell2mat(struct2cell(defaultsp));
                sptable=table('RowNames',fieldnames(defaultsp));
                sptable.grid=spvalues;
                sptable.traj=spvalues;
                app.SolverParTable.Data=sptable;
                
                defaultop=clODEfeatures.defaultObserverParams(); %app.clo_g.op
                optable=table('RowNames',fieldnames(defaultop));
                optable.value=cell2mat(struct2cell(defaultop));
                app.gridObserverParTable.Data=optable;
                app.gridObserverParTable.RowName=fieldnames(defaultop);
                
                app.featureDropDown.Items=app.clo_g.featureNames();
                app.fscale=1;
                
            else %to avoid reverting to default sp and op:
                app.clo_g.setNewProblem(app.prob);
                app.clo_t.setNewProblem(app.prob);
                
                %remove grid.name to allow new grid to be set below
                app.listenerGrid.Enabled=0;
                app.grid.name={'';''};
                app.listenerGrid.Enabled=1;
                app.clo_g.F=[]; %to allow new nVar
                app.clo_g.Xf=[]; %to allow new nVar
                app.XF=[]; %to allow new nVar
            end
            app.buildCL('grid');
            app.clo_g.getNFeatures(); %update features (nVar may affect)
            app.featureDropDown.Items=app.clo_g.featureNames();
                
            app.buildCL('traj');
            
            %solver opts from ODE file
            dt=app.prob.opt.dt;
            %dtmax, abstol/reltol, nout, maxstore?
            app.SolverParTable.Data{'dt',:}=dt;
            
            %tspan
            t0=app.prob.opt.t0;
            tf=app.prob.opt.total;
            app.clo_g.tspan=[t0,tf];
            app.clo_t.tspan=[t0,tf];
            
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
            app.gridvars=newGridvars; %triggers postSet listener

            %set up selectable gridTable name
            app.gridTable.ColumnFormat={app.gridvars.name',[],[],[],[]};
            
            % specify Z for +/- incrementing
            app.z=app.gridvars(3,:);
            app.dz=(app.z.ub-app.z.lb)/10;
            app.gridZDropDown.Items=app.gridvars.name;
            app.gridZDropDown.Value=app.z.name;
            app.dzEditField.Value=app.dz;
            
            % trajectory variables
            newTrajvars=table;
            newTrajvars.name=[app.prob.varNames(:);app.prob.auxNames(:)];
            newTrajvars.type=categorical([ones(app.prob.nVar,1);2*ones(app.prob.nAux,1)],[1,2],{'var','aux'});
            newTrajvars.ix=[(1:app.prob.nVar)';(1:app.prob.nAux)'];
            newTrajvars.Properties.RowNames=newTrajvars.name;
            app.trajvars=newTrajvars;
            
            %trajectory plot variables
            app.trajXDropDown.Items=[{'t'};app.trajvars.name];
            app.trajYDropDown.Items=app.trajvars.name;
            app.trajZDropDown.Items=[{'none'};app.trajvars.name];
            
            xlabel(app.axTrajP0,app.trajXDropDown.Value);
            ylabel(app.axTrajP0,app.trajYDropDown.Value);
            xlabel(app.axTrajClick(end),app.trajXDropDown.Value);
            ylabel(app.axTrajClick(end),app.trajYDropDown.Value);
            
            %make the actual grid data for simulation
            app.makeGridData();
            app.makeTrajData(); %p0: no clickCoords
            
            app.clo_g.initialize();
%             app.clo_t.initialize();
            
            figure(app.figControl)
        end
        
    end
    
    
    methods (Access = public)
        %constructor
        function app=gridtool(odefile)%, precision, selectedDevice, stepper, observer
            
            app.createUIComponents();
            app.attachListeners();
            
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
                app.processNewODEfile(odefile)
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
            app.figControl = figure(111);
            app.figControl.Visible = 'off';
            app.figControl.Position = [50 450 300 650];
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
            
            % Create dzEditFieldLabel
            app.dzEditFieldLabel = uilabel(app.GridTab);
            app.dzEditFieldLabel.Position = [85 510 15 20];
            app.dzEditFieldLabel.Text = 'dz';
            app.dzEditFieldLabel.HorizontalAlignment = 'right';
            
            % Create dzEditField
            app.dzEditField = uieditfield(app.GridTab, 'numeric');
            app.dzEditField.ValueChangedFcn = @app.dzEditFieldValueChanged;
            app.dzEditField.Position = [105 510 50 20];
            app.dzEditField.Value=0;
            app.dzEditField.ValueDisplayFormat='%.2g';
            
            % Create featureDropDownLabel
            app.featureDropDownLabel = uilabel(app.GridTab);
            app.featureDropDownLabel.Position = [200 510 85 20];
            app.featureDropDownLabel.Text = 'feature';
            app.featureDropDownLabel.HorizontalAlignment='center';
            
            % Create featureDropDown
            app.featureDropDown = uidropdown(app.GridTab);
            app.featureDropDown.ValueChangedFcn = @app.featureDropDownValueChanged;
            app.featureDropDown.Position = [160 490 120 20];
            
            % Create fscaleEditFieldLabel
            app.fscaleEditFieldLabel = uilabel(app.GridTab);
            app.fscaleEditFieldLabel.Position = [195 465 40 20];
            app.fscaleEditFieldLabel.Text = 'fscale';
            
            % Create fscaleEditField
            app.fscaleEditField = uieditfield(app.GridTab, 'numeric');
            app.fscaleEditField.ValueChangedFcn = @app.fscaleEditFieldValueChanged;
            app.fscaleEditField.Position = [235 465 50 20];
            app.fscaleEditField.Value=1;
            
            
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
            app.parTable.Position = [5 240 280 160];
            
            
            % Create icLabel
            app.icLabel = uilabel(app.GridTab);
            app.icLabel.Position = [5 210 91 23];
            app.icLabel.Text = 'Initial conditions';
            
            % Create icDefaultButton
            app.icDefaultButton = uibutton(app.GridTab, 'push');
            app.icDefaultButton.ButtonPushedFcn = @app.icDefaultButtonPushed;
            app.icDefaultButton.Position = [233 213 52 20];
            app.icDefaultButton.Text = 'default';
            
            % Create icTable
            app.icTable = uitable(app.GridTab);
            app.icTable.RowName={};
            app.icTable.ColumnName = {'name', 'val', 'lb', 'ub'};
            app.icTable.ColumnWidth = '1x';
            app.icTable.ColumnSortable = [true false false false];
            app.icTable.ColumnEditable = [true true true true];
            app.icTable.CellEditCallback = @app.icTableCellEdit;
            app.icTable.Position = [5 70 280 140];
            
            
            %             % Create icRandomButton
            %             app.icRandomButton = uibutton(app.GridTab, 'push');
            %             app.icRandomButton.ButtonPushedFcn = @app.icRandomButtonPushed;
            %             app.icRandomButton.Position = [229 172 52 22];
            %             app.icRandomButton.Text = 'random';
            
            
            %             app.TrajTab = uitab(app.TabGroup);
            %             app.TrajTab.Title = 'Trajectories';
            
            % Create trajectoriesLabel
            app.trajectoriesLabel = uilabel(app.GridTab);
            app.trajectoriesLabel.Position = [15 34 68 23];
            app.trajectoriesLabel.Text = 'Trajectories';
            
            % Create trajXDropDownLabel
            app.trajXDropDownLabel = uilabel(app.GridTab);
            app.trajXDropDownLabel.HorizontalAlignment = 'right';
            app.trajXDropDownLabel.Position = [15 7 10 22];
            app.trajXDropDownLabel.Text = 'x';
            
            % Create trajXDropDown
            app.trajXDropDown = uidropdown(app.GridTab);
            app.trajXDropDown.ValueChangedFcn = @app.trajXDropDownValueChanged;
            app.trajXDropDown.Position = [32 7 52 21];
            
            % Create trajYDropDownLabel
            app.trajYDropDownLabel = uilabel(app.GridTab);
            app.trajYDropDownLabel.HorizontalAlignment = 'right';
            app.trajYDropDownLabel.Position = [83 7 10 22];
            app.trajYDropDownLabel.Text = 'y';
            
            % Create trajYDropDown
            app.trajYDropDown = uidropdown(app.GridTab);
            app.trajYDropDown.ValueChangedFcn = @app.trajYDropDownValueChanged;
            app.trajYDropDown.Position = [100 7 52 21];
            
            % Create trajZDropDownLabel
            app.trajZDropDownLabel = uilabel(app.GridTab);
            app.trajZDropDownLabel.HorizontalAlignment = 'right';
            app.trajZDropDownLabel.Position = [152 7 10 22];
            app.trajZDropDownLabel.Text = 'z';
            
            % Create trajZDropDown
            app.trajZDropDown = uidropdown(app.GridTab);
            app.trajZDropDown.ValueChangedFcn = @app.trajZDropDownValueChanged;
            app.trajZDropDown.Position = [169 7 52 21];
            
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
            app.figGrid = figure(112);
            app.figGrid.Visible = 'off';
            app.figGrid.Position = [375 450 600 600];
            app.figGrid.Name = 'grid';
            app.figGrid.NumberTitle='off';
            
            % grid axis
            t=tiledlayout(app.figGrid,1,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axGrid = nexttile(t);
            
            app.imGrid=imagesc(app.axGrid,0, 0, 0);
            app.imGrid.HitTest='off';  %not clickable
            
            app.axGrid.YDir='normal';
            axis(app.axGrid,'tight');
            app.gridCBar=colorbar('northoutside');
            
            hold on
            app.markerP0=plot(nan,nan,'ko');
            app.markerP0.ButtonDownFcn=@app.clickP0;  %p0 marker is clickable
            app.markerPquery=plot(nan,nan,'ko');
            app.markerPquery.HitTest='off'; %not clickable
            hold off
            
            xlabel(app.axGrid,'x')
            ylabel(app.axGrid,'y')
            
            app.axGrid.ButtonDownFcn=@app.clickTraj; %axis is clickable
            app.figGrid.KeyPressFcn = @app.gridKeyPress;
            app.figGrid.Visible = 'on';
        end
        
        function createTrajP0Fig(app)
            % Create figTraj and hide until all components are created
            app.figTrajP0 = figure(113);
            app.figTrajP0.Visible = 'off';
            app.figTrajP0.Position = [375 250 600 200];
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
            xlabel(app.axyyTrajP0,'t')
            
            app.line3DTrajP0=plot3(app.ax3DTrajP0,nan,nan,nan);
            app.ax3DTrajP0.Visible='off';
            app.figTrajP0.Visible = 'on';
        end
        
        function createTrajClickFig(app)
            % Create figTraj and hide until all components are created
            app.figTrajClick = figure(114);
            app.figTrajClick.Visible = 'off';
            app.figTrajClick.Position = [1000 50 900 950];
            app.figTrajClick.Name = 'click trajectories';
            app.figTrajClick.NumberTitle='off';
            
            app.tilesTrajClick=tiledlayout(app.figTrajClick,app.nClick,1);
            app.tilesTrajClick.TileSpacing = 'compact';
            app.tilesTrajClick.Padding = 'compact';
            for i=1:app.nClick
                app.axTrajClick(i) = nexttile(app.tilesTrajClick);
                app.axTrajClick(i).XTickLabel=[];
            end
            app.axTrajClick(end).XTickLabelMode='auto';
            xlabel(app.axTrajClick(3),'t')
            ylabel(app.axTrajClick(3),'x')
            
            %             app.figTrajClick.Visible = 'on';
        end
    end
    
end