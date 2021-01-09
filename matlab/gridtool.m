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
        
        prob %read-only copy of the parseODE output (default values)
        %.par, parNames
        %.var, varNames
        
        devices
        grid_device
        clo_g %clODE grid feature solver
        traj_device
        clo_t %clODE trajectory solver
        
        nGrid=[32;32;10] %nPts in x, y, z dims: x,y is the 2D integration plane
        grid %array of struct encapsulating the 2-D grid properties
        %.name = valid parName/varName
        %.N = grid point values
        icNames %varNames+'0'
        gridvar2ix %mapping of par/var to their array indices
        gridvartype %mapping from par/var name to type - par or var.
        
        tspan
        nClick=3
        
        %mechanism to cache F(:,:,z) values - for rapidly scanning z post compute
        
    end
    
    % Properties that correspond to UI components
    properties (Access = public)
        
        %UI state
        gridclInitialized=false %flag to indicate if ready for integration
        trajclInitialized=false
        solutionIsCurrent=false %flag to indicate need for new solution
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %grid figure
        figGrid
        axGrid %axis handle
        imGrid %image handle
        markerP0
        markerPquery %markers specifying (x,y) for n query trajectories

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figTrajP0 %figure with p0 trajectory
        axTrajP0
        trajP0 %trajectory struct for p0 - support dual y-axis
        %.x,y, ...
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %click trajectory figure
        figTrajClick %figure with n trajectories
        axTrajClick %array of axes
        trajClick %array of trajectory struct for pQuery
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %control figure
        figControl              matlab.ui.Figure
        TabGroup                matlab.ui.container.TabGroup
        
        NumericsTab             matlab.ui.container.Tab
        currentFileLabel        matlab.ui.control.Label
        SelectODEfileButton     matlab.ui.control.Button
        gridDeviceDropDownLabel   matlab.ui.control.Label
        gridDeviceDropDown      matlab.ui.control.DropDown
        gridPrecisionCheckBox   matlab.ui.control.CheckBox
        gridStepperDropDownLabel  matlab.ui.control.Label
        gridStepperDropDown     matlab.ui.control.DropDown
        gridSolverParLabel      matlab.ui.control.Label
        gridSolverParTable      matlab.ui.control.Table
        gridObserverDropDownLabel   matlab.ui.control.Label
        gridObserverDropDown    matlab.ui.control.DropDown
        gridObserverParTable    matlab.ui.control.Table
        gridObserverParLabel    matlab.ui.control.Label
        
        trajSolverParTable      matlab.ui.control.Table
        trajSolverParLabel      matlab.ui.control.Label
        trajDeviceDropDownLabel   matlab.ui.control.Label
        trajStepperDropDownLabel  matlab.ui.control.Label
        trajStepperDropDown     matlab.ui.control.DropDown
        trajPrecisionCheckBox   matlab.ui.control.CheckBox
        trajDeviceDropDown      matlab.ui.control.DropDown
        
        GridTab                 matlab.ui.container.Tab
        featureDropDownLabel    matlab.ui.control.Label
        featureDropDown         matlab.ui.control.DropDown
        fscaleEditFieldLabel    matlab.ui.control.Label
        fscaleEditField         matlab.ui.control.NumericEditField
        t0EditFieldLabel        matlab.ui.control.Label
        t0EditField             matlab.ui.control.NumericEditField
        tfEditFieldLabel        matlab.ui.control.Label
        tfEditField             matlab.ui.control.NumericEditField
        gridTable               matlab.ui.control.Table
        parLabel                matlab.ui.control.Label
        parTable                matlab.ui.control.Table
        parDefaultButton        matlab.ui.control.Button
        icLabel                 matlab.ui.control.Label
        icTable                 matlab.ui.control.Table
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
        function gridDeviceDropDownValueChanged(app, event)
            value = app.gridDeviceDropDown.Value;
            app.clo_g.setOpenCL(value);
            app.gridclInitialized=false;
        end

        % Value changed function: trajDeviceDropDown
        function trajDeviceDropDownValueChanged(app, event)
            value = app.trajDeviceDropDown.Value;
            app.clo_t.setOpenCL(value); 
            app.trajclInitialized=false;           
        end

        % Value changed function: gridPrecisionCheckBox
        function gridPrecisionCheckBoxValueChanged(app, event)
            useDouble = app.gridPrecisionCheckBox.Value;
            precision='single';
            if useDouble, precision='double'; end
            app.clo_g.setPrecision(precision);
            app.gridclInitialized=false;
        end

        % Value changed function: trajPrecisionCheckBox
        function trajPrecisionCheckBoxValueChanged(app, event)
            useDouble = app.trajPrecisionCheckBox.Value;
            precision='single';
            if useDouble, precision='double'; end
            app.clo_t.setPrecision(precision);
            app.trajclInitialized=false;
        end

        % Value changed function: gridStepperDropDown
        function gridStepperDropDownValueChanged(app, event)
            value = app.gridStepperDropDown.Value;
            
        end

        % Value changed function: trajStepperDropDown
        function trajStepperDropDownValueChanged(app, event)
            value = app.trajStepperDropDown.Value;
            
        end

        % Value changed function: gridObserverDropDown
        function gridObserverDropDownValueChanged(app, event)
            value = app.gridObserverDropDown.Value;
            
        end

        % Cell edit callback: gridSolverParTable
        function gridSolverParTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Cell edit callback: gridObserverParTable
        function gridObserverParTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Cell edit callback: trajSolverParTable
        function trajSolverParTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Value changed function: featureDropDown
        function featureDropDownValueChanged(app, event)
            value = app.featureDropDown.Value;
            
        end

        % Value changed function: fscaleEditField
        function fscaleEditFieldValueChanged(app, event)
            value = app.fscaleEditField.Value;
            
        end

        % Value changed function: t0EditField
        function t0EditFieldValueChanged(app, event)
            value = app.t0EditField.Value;
            
        end

        % Value changed function: tfEditField
        function tfEditFieldValueChanged(app, event)
            value = app.tfEditField.Value;
            
        end

        % Cell edit callback: gridTable
        function gridTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Cell edit callback: parTable
        function parTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Cell edit callback: icTable
        function icTableCellEdit(app, event)
            indices = event.Indices;
            newData = event.NewData;
            
        end

        % Button pushed function: parDefaultButton
        function parDefaultButtonPushed(app, event)
            
        end

        % Button pushed function: icDefaultButton
        function icDefaultButtonPushed(app, event)
            
        end

        % Button pushed function: icRandomButton
        function icRandomButtonPushed(app, event)
            
        end

        % Value changed function: xDropDown
        function xDropDownValueChanged(app, event)
            value = app.xDropDown.Value;
            
        end

        % Value changed function: yDropDown
        function yDropDownValueChanged(app, event)
            value = app.yDropDown.Value;
            
        end

        % Value changed function: zDropDown
        function zDropDownValueChanged(app, event)
            value = app.zDropDown.Value;
            
        end

        % Value changed function: clicksEditField
        function clicksEditFieldValueChanged(app, event)
            value = app.clicksEditField.Value;
            
        end

        % Value changed function: linkAxesButton
        function linkAxesButtonValueChanged(app, event)
            value = app.linkAxesButton.Value;
            
        end

        % Value changed function: threeDButton
        function threeDButtonValueChanged(app, event)
            value = app.threeDButton.Value;
            
        end
    end
    
    % UI interaction callbacks
    methods (Access = private)
    end
    
    % helper functions
    methods (Access = private)
        
        function processNewODEfile(app, odefile, firstTime)
            
            if ~exist('firstTime','var')
                firstTime=false;
            end
            
            app.odefile=odefile;
            [~,app.prob]=ode2cl(odefile);
            
            if app.prob.nPar+app.prob.nVar<2
                error('Not enough parameters and variables for a 2D grid!')
            end
            app.icNames=strcat(app.prob.varNames(:),'0');
            app.gridvar2ix=containers.Map([app.prob.parNames(:);...
                app.icNames(:)],...% app.prob.varNames(:)],...
                [(1:app.prob.nPar)';(1:app.prob.nVar)']);
            
            %If first time, set up clODE objects using their default parameters
            if firstTime
                app.clo_g=clODEfeatures(app.prob, [], app.grid_device);
                app.gridObserverDropDown.Items=app.clo_g.observerNames;
                app.gridObserverDropDown.Value=app.clo_g.observer;
                
                app.clo_t=clODEtrajectory(app.prob, [], app.traj_device);
            
                app.gridStepperDropDown.Items=app.clo_g.getAvailableSteppers();
                app.gridStepperDropDown.Value=app.clo_g.stepper;
                app.trajStepperDropDown.Items=app.clo_t.getAvailableSteppers();
                app.trajStepperDropDown.Value=app.clo_t.stepper;
            
                sptable=table;
                sptable.name=fieldnames(app.clo_g.sp);
                sptable.value=cell2mat(struct2cell(app.clo_g.sp));
                app.gridSolverParTable.Data=sptable;
                app.trajSolverParTable.Data=sptable;
            
                optable=table;
                optable.name=fieldnames(app.clo_g.op);
                optable.value=cell2mat(struct2cell(app.clo_g.op));
                app.gridObserverParTable.Data=optable;
                
                app.featureDropDown.Items=app.clo_g.featureNames();
                
                app.grid=table('Size',[3,3],'VariableTypes',{'categorical','double','double'},...
                    'VariableNames',{'name','N','del'},'RowNames',{'x';'y';'z'}); 
                app.grid.N=app.nGrid;
            else
                app.clo_g.setNewProblem(app.prob);
                app.clo_t.setNewProblem(app.prob);
            end
            app.gridclInitialized=false;
            app.trajclInitialized=false;
            
            %populate UI components
            app.currentFileLabel.Text=app.prob.name+".ode";
            
            app.t0EditField.Value=app.prob.opt.t0;
            app.tfEditField.Value=app.prob.opt.total;
            
            app.parTable.Data=struct2table(app.prob.par);
            
            fullvartable=struct2table(app.prob.var);
            app.icTable.Data=fullvartable(:,1:4);
            
            %set up initial default grid
            if app.prob.nPar<1
                gridnames{1}=app.prob.var(1).name;
                gridnames{2}=app.prob.var(2).name;
                gridnames{3}=app.prob.var(3).name;
            elseif app.prob.nPar<2
                gridnames{1}=app.prob.par(1).name;
                gridnames{2}=app.prob.var(1).name;
                gridnames{3}=app.prob.var(2).name;
            elseif app.prob.nPar<3
                gridnames{1}=app.prob.par(1).name;
                gridnames{2}=app.prob.par(2).name;
                gridnames{3}=app.prob.var(1).name;
            else
                gridnames{1}=app.prob.par(1).name;
                gridnames{2}=app.prob.par(2).name;
                gridnames{3}=app.prob.par(3).name;
            end
            
            app.axGrid.XAxis.Label.String=gridnames{1};
            app.axGrid.YAxis.Label.String=gridnames{2};
            
            allgridvarnames = [app.prob.parNames';app.icNames];
            app.grid.name=categorical(gridnames(:),allgridvarnames);
            
            
            app.gridTable.Data=app.grid;
            
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
            app.figControl.Position = [50 400 300 600];
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
            app.gridSolverParTable.ColumnName = {'name'; 'value'};
            app.gridSolverParTable.ColumnWidth = {'auto', 'fit'};
            app.gridSolverParTable.RowName = {};
            app.gridSolverParTable.ColumnEditable = [false true];
            app.gridSolverParTable.CellEditCallback = @app.gridSolverParTableCellEdit;
            app.gridSolverParTable.Position = [13 355 150 145];
            

            % Create gridObserverParLabel
            app.gridObserverParLabel = uilabel(app.NumericsTab);
            app.gridObserverParLabel.Position = [13 330 119 22];
            app.gridObserverParLabel.Text = 'Observer parameters';

            % Create gridObserverParTable
            app.gridObserverParTable = uitable(app.NumericsTab);
            app.gridObserverParTable.ColumnName = {'name'; 'value'};
            app.gridObserverParTable.ColumnWidth = {'auto', 'fit'};
            app.gridObserverParTable.RowName = {};
            app.gridObserverParTable.ColumnEditable = [false true];
            app.gridObserverParTable.CellEditCallback = @app.gridObserverParTableCellEdit;
            app.gridObserverParTable.Position = [13 186 150 145];

            % Create gridObserverDropDownLabel
            app.gridObserverDropDownLabel = uilabel(app.NumericsTab);
            app.gridObserverDropDownLabel.HorizontalAlignment = 'center';
            app.gridObserverDropDownLabel.Position = [172 330 52 22];
            app.gridObserverDropDownLabel.Text = 'observer';

            % Create gridObserverDropDown
            app.gridObserverDropDown = uidropdown(app.NumericsTab);
            app.gridObserverDropDown.ValueChangedFcn = @app.gridObserverDropDownValueChanged;
            app.gridObserverDropDown.Position = [175 309 107 22];
            

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
            app.trajSolverParTable.ColumnName = {'name'; 'value'};
            app.trajSolverParTable.ColumnWidth = {'auto', 'fit'};
            app.trajSolverParTable.RowName = {};
            app.trajSolverParTable.ColumnEditable = [false true];
            app.trajSolverParTable.CellEditCallback = @app.trajSolverParTableCellEdit;
            app.trajSolverParTable.Position = [15 15 150 145];
            

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
            
            % Create t0EditFieldLabel
            app.t0EditFieldLabel = uilabel(app.GridTab);
            app.t0EditFieldLabel.HorizontalAlignment = 'center';
            app.t0EditFieldLabel.Position = [176 498 43 22];
            app.t0EditFieldLabel.Text = 't0';

            % Create t0EditField
            app.t0EditField = uieditfield(app.GridTab, 'numeric');
            app.t0EditField.ValueChangedFcn = @app.t0EditFieldValueChanged;
            app.t0EditField.Position = [220 500 55 18];

            % Create tfEditFieldLabel
            app.tfEditFieldLabel = uilabel(app.GridTab);
            app.tfEditFieldLabel.HorizontalAlignment = 'center';
            app.tfEditFieldLabel.Position = [176 479 43 22];
            app.tfEditFieldLabel.Text = 'tf';

            % Create tfEditField
            app.tfEditField = uieditfield(app.GridTab, 'numeric');
            app.tfEditField.ValueChangedFcn = @app.tfEditFieldValueChanged;
            app.tfEditField.Position = [220 481 55 18];

            % Create gridTable
            app.gridTable = uitable(app.GridTab);
            app.gridTable.ColumnName = {'name'; 'N'; 'dvar'};
            app.gridTable.ColumnWidth = {'auto', 'fit', 'fit'};
            app.gridTable.ColumnEditable = [true true true];
            app.gridTable.CellEditCallback = @app.gridTableCellEdit;
            app.gridTable.Position = [15 430 150 90];
            

            % Create parLabel
            app.parLabel = uilabel(app.GridTab);
            app.parLabel.Position = [15 403 69 23];
            app.parLabel.Text = 'Parameters';

            % Create parTable
            app.parTable = uitable(app.GridTab);
            app.parTable.ColumnName = {'name'; 'value'; 'lb'; 'ub'};
            app.parTable.ColumnWidth = {'auto', 'fit', 'fit', 'fit'};
            app.parTable.RowName = {};
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
            app.icTable.ColumnName = {'name'; 'value'; 'lb'; 'ub'};
            app.icTable.ColumnWidth = {'auto', 'fit', 'fit', 'fit'};
            app.icTable.RowName = {};
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
            
            app.figControl.Visible = 'on';
            
        end
            
        function createGridFig(app)
            % Create figGrid and hide until all components are created
            app.figGrid = uifigure('Visible', 'off');
            app.figGrid.Position = [375 400 600 600];
            app.figGrid.Name = 'grid';
            
            % grid axis
            t=tiledlayout(app.figGrid,1,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axGrid = nexttile(t);
            xlabel(app.axGrid,'x')
            ylabel(app.axGrid,'y')
%             app.axGrid = axes(app.figGrid);
            
            app.figGrid.Visible = 'on';
        end
        
        function createTrajP0Fig(app)
            % Create figTraj and hide until all components are created
            app.figTrajP0 = uifigure('Visible', 'off');
            app.figTrajP0.Position = [50 50 925 300];
            app.figTrajP0.Name = 'p0';
            
            t=tiledlayout(app.figTrajP0,1,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axTrajP0 = nexttile(t);
            xlabel(app.axTrajP0,'t')
            ylabel(app.axTrajP0,'x')
%             app.axTrajP0 = axes(app.figTrajP0);
            
            app.figTrajP0.Visible = 'on';
        end
        
        function createTrajClickFig(app)
            % Create figTraj and hide until all components are created
            app.figTrajClick = uifigure('Visible', 'off');
            app.figTrajClick.Position = [1000 50 900 950];
            app.figTrajClick.Name = 'click trajectories';

            t=tiledlayout(app.figTrajClick,3,1);
            t.TileSpacing = 'compact';
            t.Padding = 'compact';
            app.axTrajClick(1) = nexttile(t);
            app.axTrajClick(2) = nexttile(t);
            app.axTrajClick(3) = nexttile(t);
            xlabel(app.axTrajClick(3),'t')
            ylabel(app.axTrajClick(3),'x')

            app.figTrajClick.Visible = 'on';
        end
    end
    
end