classdef gridtool_old<handle
    
    %TODO:
    % Build panel:
    % - device selection: update clODE to have traj-specific solverpars
    % (for both grid + traj on one device); allow separate device for grid
    % vs traj
    % - select ode file, stepper, observer, precision; BUILD button
    % - allow passing in values for all above
    %
    % Parameters/ICs:
    % - save parameters AND ICs to file (ORIGIN or Current set)
    % - show/set/save p --> show/set/save ORIGIN (fits better under GRID controls)
    %
    % Grids
    % - multi-window support (e.g. different features)
    %
    % Trajectories:
    % - add function to simulate current point (ORIGIN)'s trajectory
    % - how to handle multi-click vs. ORIGIN trajectory?
    % - separate control real points vs snap2grid
    % - multi-window support
    % - hotkeys for Window Fit, Viewaxes 2D/3D, Make Create, etc...
    %
    % General cleanup
    % - subroutine to check numeric/textual inputs (called when any textbox/uitable is edited)
    % - better way to manage linking changes across uitables: notifier/listeners?
    % - define grid/traj classes OR use structs to organize data better? go full OO? 
    
    %TODO: make most of these things private, create get/set functions
    %where appropriate
    properties %(SetAccess=private, GetAccess = private)
        
        %problem data
        p       %current default parameters
        pNames
        nPar
        plb
        pub
        
        y0      %current default initial condition
        vNames
        v0Names %add "_0" to vNames
        nVar
        y0lb
        y0ub
        
        auxNames
        
        %grid - could make a class for this. getFeatures() would
        %hide method of commputation, allowing different backends (eg. clODE)
        xIx
        xispar 
        gridXname
        nx
        yIx
        yispar
        gridYname
        ny
        nPts    %nPts=nx*ny
        
        %gridded params
        x    %grid as vector 
        xx   %from meshgrid
        y
        yy
        f
        ff
        
        %clODE grid data
        Pars    %parameter matrix
        Y0      %initial condition matrix
        Yf      %final solution values
        F       %Features matrix
        
        fFunc=@(f,n) deal(f, n);
        fNames;
        
        %trajectory - could make a class for this. getTrajectories() would
        %hide method of commputation, allowing different backends (eg. clODE)
        t
        trajPars
        trajY0
        trajY0def %for 'go'
        traj
        trajDY
        trajAux
        trajYf
        trajTF
        
        %solver + observer parameters
        dt
        Times
        stepper
        observer
        obspars
        
        dtTraj
        TimesTraj
        
        nout
        
        %original problem data/options
        clFilename
        odeFilename
        ivp
        XPPopt
        gridTableDefault
        parTableDefault
        y0TableDefault
        numericsTableDefault
        obsparDefault
        obsparNames
        
        %control window with uitables & callbacks to set sim data
        controlFig
        gridTable
        parTable
        y0Table
        numericsTable
        observerTable
        
        featureSelect
        gridVariableSelect
        
        zParSelect
        zIx
        dz
        dzEdit
        zValueText
        
        plotControlPanel
        featureScaleSelect
        plottypeSelect
        cAxisPanel
        CAxisMinEdit
        CAxisMaxEdit
        CAxisFitButton
        
        trajControlPanel
        trajXSelect
        trajXnames
        trajYSelect
        trajYnames
        trajZSelect
        trajZnames
        trajTable
        trajTableDefault
            
        numericsTableTraj
        numericsTableDefaultTraj
        trajNumericsNames
        gridNumericsNames
                
        %trajectory figure. 
        trajFig    
        trajAxisMode='auto';
        nClick=3;
        xClickVal
        yClickVal
        markers='osdp^v';
        markerColor='k';
        trajSubPlotAx
        
        %gridplot stuff
        gridFig   
        gridAxis
        selectedPlotType='image';
        availablePlotTypes={'image','scatter','contourf'};
        featIx=1;
        gridPlotHandle
        showParToggle
        cmin
        cmax
        cAxisMode='auto';
        
        %UI state
        lastkeypress
        havePreviousGrid=false;
        havePreviousTraj=false;
        wasTransient=false;
        gridIsDirty=true; %indicate if plot and problem data are out of sync
        
        %clODE objects
        precision='single';
        clODEgrid
        clODEtraj
        
    end
    
    methods
        
        %constructor
        function gt=gridtool(clPlatform,clDevice,odefile,stepper,observer,newFfunc)
            
            %set defaults of inputs not provided
            if ~exist('odefile','var')||isempty(odefile)
                [name,path]=uigetfile('.ode','Select an ODE file');
                if ~ischar(name)
                    disp('File selection canceled, quitting...')
                    return
                end
                odefile=fullfile(path,name);
            end
            
            if ~exist('stepper','var')||isempty(stepper)
                stepper='rungekutta4';
            end
            gt.stepper=stepper;
            
            if ~exist('observer','var')||isempty(observer)
                observer='localmin';
            end
            gt.observer=observer;
            
            if ~exist('newFfunc','var')||isempty(newFfunc)
                newFfunc=gt.fFunc;
            end
            gt.fFunc=newFfunc;
            
            if ~exist('clPlatform','var')||isempty(clPlatform)
            end
            if ~exist('clDevice','var')||isempty(clDevice)
            end
            
            %convert odefile to cl, and extract problem data
%             [gt.clFilename, gt.p, gt.y0, gt.ivp, gt.XPPopt]=ode2cl(odefile,gt.precision);
            [~, xppdata]=xppConverter(odefile,'cl');
            
            gt.ivp=xppdata;
            gt.XPPopt=xppdata.opt;
            gt.p=xppdata.p0;
            gt.y0=xppdata.x0;
            
            gt.processParserOutput();
            
            %set up clODE
            clPlatformGrid=clPlatform(1);
            clDeviceGrid=clDevice(1);
            if length(clPlatform)==1
                clPlatformTraj=clPlatformGrid;
                clDeviceTraj=clDeviceGrid;
            elseif length(clPlatform)==2
                clPlatformTraj=clPlatform(2);
                clDeviceTraj=clDevice(2);
            end
            
            %keep a separate object for each because currently solverPars
            %are reused for both grid and trajectory. *change c++ class so
            %there are separate controls for traj/grid?? (initialize to
            %identical) would allow this GUI to build solely on one device,
            %reducing up setup time. (want to still allow separate device for traj)
            tic
            gt.clODEgrid=clODE(gt.ivp, stepper, observer, clPlatformGrid, clDeviceGrid, gt.fFunc);
            gt.clODEtraj=clODE(gt.ivp, stepper, observer, clPlatformTraj, clDeviceTraj, gt.fFunc);
            buildTime=toc
            
            [~,gt.fNames]=gt.clODEgrid.getF();
            
            gt.clODEgrid.setSolverParams(gt.dt, gt.Times);
            gt.clODEgrid.setObserverParams(gt.obspars);
            
            gt.clODEtraj.setSolverParams(gt.dt, gt.TimesTraj,gt.nout);
            gt.clODEtraj.setObserverParams(gt.obspars);
            
            gt.generateProblemData();
            
            gt.buildFigs();

            figure(gt.gridFig);
        end
        
        
    end
    
    
    %public methods
    methods
        
        
    function setPars(gt,p)

        gt.p=p(:);

        tempParTableData=get(gt.parTable,'Data');
        tempParTableData(:,2)=num2cell(p(:));
        set(gt.parTable,'Data',tempParTableData);

        set(gt.zValueText,'String',[gt.pNames{gt.zIx} '=' num2str(gt.p(gt.zIx))]);
        gt.generateProblemData();
        gt.toggleShowParsMarker();
    end
    
    
    end
    
    
    %plotting, etc
    methods (Access=private)
        
        
        function processParserOutput(gt)
            
            gt.odeFilename=gt.ivp.name;
        
            gt.pNames={gt.ivp.par(:).name};
            gt.nPar=length(gt.p);
            gt.plb=[gt.ivp.par(:).lb];
            gt.pub=[gt.ivp.par(:).ub];
            
            gt.vNames={gt.ivp.var(:).name};
            gt.nVar=length(gt.y0);
            gt.y0lb=[gt.ivp.var(:).lb];
            gt.y0ub=[gt.ivp.var(:).ub];
            gt.v0Names=strcat(gt.vNames,repmat({'_0'},1,gt.nVar));
            
            gt.auxNames={gt.ivp.aux(:).name};
            
            %expand non-ranged params to +- 50%
            notRangedP= gt.pub-gt.plb == 0;
            gt.plb(notRangedP)=gt.p(notRangedP)-0.5*abs(gt.p(notRangedP));
            gt.pub(notRangedP)=gt.p(notRangedP)+0.5*abs(gt.p(notRangedP)); 
            
            notRangedy0= gt.y0ub-gt.y0lb == 0;
            gt.y0lb(notRangedy0)=gt.y0(notRangedy0)-0.25*abs(gt.y0(notRangedy0));
            gt.y0ub(notRangedy0)=gt.y0(notRangedy0)+0.25*abs(gt.y0(notRangedy0));
            
            %if default p value is zero, make range nontrivial
            gt.pub(gt.pub==gt.plb)=gt.pub(gt.pub==gt.plb)+1;
            gt.y0ub(gt.y0ub==gt.y0lb)=gt.y0ub(gt.y0ub==gt.y0lb)+1;
            
            %default grid on first two parameters
            gt.xIx=1;
            gt.gridXname=gt.pNames{gt.xIx};
            gt.xispar=true;
            gt.nx=32;
            
            gt.yIx=2;
            gt.gridYname=gt.pNames{gt.yIx};
            gt.yispar=true;
            gt.ny=32;
            
            gt.nPts=gt.nx*gt.ny;
            
            % z-par select
            if gt.nPar>2
                gt.zIx=3;
            else
                gt.zIx=2; %???
            end
            gt.dz=(gt.pub(gt.zIx)-gt.plb(gt.zIx))/10;
            
            
            %set the default solver+observer params so sim is ready to go
            gt.dt=gt.XPPopt.dt;
            gt.dtTraj=gt.XPPopt.dt;
            
            if gt.XPPopt.t0==gt.XPPopt.trans
                gt.XPPopt.trans=(gt.XPPopt.total-gt.XPPopt.t0)/2;
            end
            gt.Times=[gt.XPPopt.t0,gt.XPPopt.trans,gt.XPPopt.total];
            gt.TimesTraj=gt.Times;
            
            gt.nout=1;

            %feature observer parameters (not all observers use all of these...)
            gt.obspars=observerset();
            
            %Use defaults, or get from XPP file??
%             gt.obspars.minYamp=gt.XPPopt.minYamp;
%             gt.obspars.minDYamp=gt.XPPopt.minDYamp;
%             gt.obspars.fractionYup=gt.XPPopt.fractionYup;
%             gt.obspars.fractionYdown=gt.XPPopt.fractionYdown;
%             gt.obspars.fractionDYup=gt.XPPopt.fractionDYup;
%             gt.obspars.fractionDYdown=gt.XPPopt.fractionDYdown;
%             gt.obspars.minIMI=gt.XPPopt.minIMI;
%             gt.obspars.normTol=gt.XPPopt.normTol;
            
            gt.obsparDefault=struct2cell(gt.obspars);
            gt.obsparNames=fieldnames(gt.obspars);
            
            switch lower(gt.observer)
                case {'plateau','plateaunoise'}
                    obsIxList=2:9;
                case {'localmin','localmax'}
                    obsIxList=[2:5,10];
            end
            gt.obsparDefault=gt.obsparDefault(obsIxList);
            gt.obsparNames=gt.obsparNames(obsIxList);
            
            dvNames=strcat('d',gt.vNames(:));
            
            gt.trajXnames=['t';gt.vNames(:);dvNames;gt.auxNames(:)];
            gt.trajYnames=[gt.vNames(:);dvNames;gt.auxNames(:)];
            gt.trajZnames=['none';gt.vNames(:);dvNames;gt.auxNames(:)];
            
            %store the default values for UI elements
            
            gt.gridTableDefault={gt.pNames{gt.xIx},gt.plb(gt.xIx),gt.pub(gt.xIx),gt.nx;... 
                gt.pNames{gt.yIx},gt.plb(gt.yIx),gt.pub(gt.yIx),gt.ny};
            gt.parTableDefault=[gt.pNames(:), num2cell([gt.p(:),gt.plb(:),gt.pub(:)])];
            gt.y0TableDefault=[gt.v0Names(:),num2cell([gt.y0(:),gt.y0lb(:),gt.y0ub(:)])];
            gt.numericsTableDefault=num2cell([gt.dt;gt.Times(:)]);
            gt.numericsTableDefaultTraj=num2cell([gt.nClick;gt.dtTraj;gt.TimesTraj([1,3])';gt.nout]); 
            
            gt.gridNumericsNames={'dt','t0','trans','total'};
            gt.trajNumericsNames={'nClick','dt','t0','total','nOut'};
            
        end
        
        
        
        function buildFigs(gt)
            
            gt.buildGridFig()
            gt.buildTrajFig()
            gt.buildControlFig()
            
        end
        
        %%%%main grid window
        function buildGridFig(gt)
            gt.gridFig=figure('NumberTitle','off',...
                'toolbar','figure',...
                'Name',['gridtool: ' gt.ivp.name],...
                'Units','normalized',...
                'Position',[.01 .35 .4 .55],...
                'WindowKeyPressFcn',@gt.processKeyPressGrid);

            gt.gridAxis = axes('Parent',gt.gridFig,...           
                'buttondownfcn',@gt.clickTrajectory);
            
            axis([gt.plb(gt.xIx) gt.pub(gt.xIx) gt.plb(gt.yIx) gt.pub(gt.yIx)]);
        end
        
        %%%%trajectory window - not visible until needed
        function buildTrajFig(gt)
            gt.trajFig=figure('NumberTitle','off',...
                'toolbar','figure',...
                'Name','gridtool: trajectories',...
                'visible','off',...
                'WindowKeyPressFcn',@gt.processKeyPressTrajectories);
        end
          
        %%%%control window - closing this one closes the whole UI
        function buildControlFig(gt)
            gt.controlFig=figure('NumberTitle','off',...
                'toolbar','figure',...
                'Name','gridtool control panel',...
                'Units','normalized',...
                'Position',[.72 .05 .275 .85],...
                'WindowKeyPressFcn',@gt.processKeyPressControl,...
                'CloseRequestFcn',@gt.my_closereq);
            
            %grid controls
            uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','Grid Controls',...
                'Units','normalized',...
                'Position',[.05 .9 .8 .02]);
            
            % Feature selection popup
            gt.featureSelect=uicontrol(gt.controlFig,...       
                'Style','popup','units','normalized',...
                'Position',[0.4,0.9,0.39,0.025],...'Position',[20,30,160,20],...
                'String',gt.fNames,...
                'Value',gt.featIx,...
                'Callback',@gt.changeFeature);
            
            % Variable selection popup
            gt.gridVariableSelect=uicontrol(gt.controlFig,...       
                'Style','popup','units','normalized',...
                'Position',[0.8,0.9,0.15,0.025],...'Position',[20,30,160,20],...
                'String',gt.vNames,...
                'Value',gt.clODEgrid.Obs.varIx,...
                'Callback',@gt.changeGridVariable);
            
            
            gt.gridTable=uitable('Parent',gt.controlFig,...
                'Units','normalized',...
                'Position',[0.05,0.8,0.9,0.1],...
                'RowName',{'x','y'}, ...
                'ColumnName', {'name','min','max','nPts'},...
                'ColumnWidth',{125,'auto','auto','auto'},...
                'ColumnFormat',{[gt.pNames,gt.v0Names],'numeric','numeric','numeric'},...
                'ColumnEditable', [true true true true],...
                'Data',gt.gridTableDefault,...
                'CellEditCallback',@gt.editGridControls);
%             gt.gridTable.Position(3) = gt.gridTable.Extent(3);

            %z: easy stepping up and down in a third dimension
            % set a deltaZ and hook up the +/- keys to step in z
            uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','z: ',...
                'Units','normalized',...
                'Position',[.075 .775 .1 .02]);
            
            gt.zParSelect=uicontrol(gt.controlFig,...       
                'Style','popup','units','normalized',...
                'Position',[0.115,0.775,0.26,0.025],...'Position',[20,30,160,20],...
                'String',gt.pNames,...
                'Value',gt.zIx,...
                'Callback',@gt.changeZpar);
            
            uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','dz: ',...
                'Units','normalized',...
                'Position',[.38 .775 .1 .02]);
            
            gt.dzEdit=uicontrol(gt.controlFig,...
                'style','edit','units','normalized',...
                'position',[0.43,0.775,0.1,0.025],...
                'String',num2str(gt.dz),...
                'Callback',@gt.editDZ);
            
            gt.zValueText=uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String',[gt.pNames{gt.zIx} '=' num2str(gt.p(gt.zIx))],...
                'Units','normalized',...
                'Position',[.55 .775 .3 .02]);
            

            %par controls
            uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','Parameters',...
                'Units','normalized',...
                'Position',[.05 .725 .8 .02]);
            
            gt.showParToggle=uicontrol('Parent',gt.controlFig,...
                'Style','togglebutton',...
                'String','show p',...
                'Units','normalized',...
                'Position',[.55 .725 .1 .025],...
                'Callback',@gt.toggleShowParsMarker);
            
            uicontrol('Parent',gt.controlFig,...
                'Style','pushbutton',...
                'String','set p',...
                'Units','normalized',...
                'Position',[.65 .725 .1 .025],...
                'Callback',@gt.clickSetPars);
            
            uicontrol('Parent',gt.controlFig,...
                'Style','pushbutton',...
                'String','save p',...
                'Units','normalized',...
                'Position',[.75 .725 .1 .025],...
                'Callback',@gt.savePars);
            
            uicontrol('Parent',gt.controlFig,...
                'Style','pushbutton',...
                'String','reset',...
                'Units','normalized',...
                'Position',[.85 .725 .1 .025],...
                'Callback',@gt.resetPars);
            
            
            gt.parTable=uitable('Parent',gt.controlFig,...
                'Units','normalized',...
                'Position',[0.05,0.575,0.9,0.15],...
                'RowName',{}, ...
                'ColumnName', {'name','value','min','max'},...
                'ColumnWidth',{150,'auto','auto','auto'},...
                'ColumnFormat',{'char','numeric','numeric','numeric'},...
                'ColumnEditable', [true true true true],...
                'Data',gt.parTableDefault,...
                'CellEditCallback',@gt.editParams);
%             gt.parTable.Position(3) = gt.parTable.Extent(3);
            
            
            %initial condition controls
            uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','Initial Conditions',...
                'Units','normalized',...
                'Position',[.05 .525 .8 .02]);
            
            uicontrol('Parent',gt.controlFig,...
                'Style','pushbutton',...
                'String','random',...
                'Units','normalized',...
                'Position',[.7 .525 .15 .025],...
                'Callback',@gt.randomizeY0);
            
            uicontrol('Parent',gt.controlFig,...
                'Style','pushbutton',...
                'String','reset',...
                'Units','normalized',...
                'Position',[.85 .525 .1 .025],...
                'Callback',@gt.resetY0);
            
            
            gt.y0Table=uitable('Parent',gt.controlFig,...
                'Units','normalized',...
                'Position',[0.05,0.4,0.9,0.125],...
                'RowName',{}, ...
                'ColumnName', {'name','value','min','max'},...
                'ColumnWidth',{150,'auto','auto','auto'},...
                'ColumnFormat',{'char','numeric','numeric','numeric'},...
                'ColumnEditable', [true true true true],...
                'Data',gt.y0TableDefault,...
                'CellEditCallback',@gt.editInitialData);
%             gt.y0Table.Position(3) = gt.y0Table.Extent(3);
            
            
            %Solver parameter controls
            uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','nUmerics',...
                'Units','normalized',...
                'Position',[.05 .35 .35 .02]);
            
            
            gt.numericsTable=uitable('Parent',gt.controlFig,...
                'Units','normalized',...
                'Position',[0.05,0.25,0.35,0.1],...
                'ColumnName', {},...
                'ColumnWidth','auto',...
                'ColumnFormat',{'numeric'},...
                'ColumnEditable', true,...
                'Data',gt.numericsTableDefault,...
                'RowName', gt.gridNumericsNames, ...
                'CellEditCallback',@gt.editSolverParams);
%             gt.numericsTable.Position(3) = gt.numericsTable.Extent(3);
            
            %Observer parameter controls
            uicontrol('Parent',gt.controlFig,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','Observer Parameters',...
                'Units','normalized',...
                'Position',[.45 .35 .45 .02]);
            
            gt.observerTable=uitable('Parent',gt.controlFig,...
                'Units','normalized',...
                'Position',[0.45,0.25,0.5,0.1],...
                'ColumnName', {},...
                'RowName', gt.obsparNames,...
                'ColumnWidth',{'auto'},...
                'ColumnFormat',{'numeric'},...
                'ColumnEditable', true,...
                'Data',gt.obsparDefault,...%'RowName', obsparNames, ...
                'CellEditCallback',@gt.editObserverParams);
%             gt.observerTable.Position(3) = gt.observerTable.Extent(3);



            gt.plotControlPanel=uipanel(gt.controlFig,...
                'Title','Grid Plot Controls','units','normalized',...
                'position',[0.01,0.01,0.49,0.22]);
            
            gt.plottypeSelect=uicontrol(gt.plotControlPanel,...          % Select the plot type
                'Style','popup','units','normalized',...
                'Position',[0.525,0.85,0.425,0.1],...'Position',[20,30,160,20],...
                'String',gt.availablePlotTypes,...
                'Value',find(strcmp(gt.availablePlotTypes,gt.selectedPlotType)),...
                'Callback',@gt.changePlotType);     
            
            %Color limits panel
            gt.cAxisPanel=uipanel(gt.plotControlPanel,...
                'Title','Color Limits','units','normalized',...
                'position',[0.05,0.6,0.9,0.225]);
            
            gt.CAxisMinEdit=uicontrol(gt.cAxisPanel,...
                'style','edit','units','normalized',...
                'position',[0.05,0.05,0.3,0.9],...
                'String','',...
                'Callback',@gt.editCmin);
            
            gt.CAxisMaxEdit=uicontrol(gt.cAxisPanel,...
                'style','edit','units','normalized',...
                'position',[0.35,0.05,0.3,0.9],...
                'String','',...
                'Callback',@gt.editCmax);
            
            gt.CAxisFitButton=uicontrol(gt.cAxisPanel,...
                'style','pushbutton','units','normalized',...
                'position',[0.65,0.025,0.3,0.95],...
                'String','Fit',...
                'Callback',@gt.fitCAxis);
            
            
            %trajectory controls
            
            gt.trajControlPanel=uipanel(gt.controlFig,...
                'Title','Trajectory Controls','units','normalized',...
                'position',[0.5,0.01,0.49,0.22]);
            
            %variable selection table
            
            uicontrol(gt.trajControlPanel,...
                'style','pushbutton','units','normalized',...
                'position',[0.07,0.89,0.25,0.1],...
                'String','Fit',...
                'Callback',@gt.fitTrajAxis);
            
            
            uicontrol('Parent',gt.trajControlPanel,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','x',...
                'Units','normalized',...
                'Position',[.01 .78 .1 .1]);
            
            gt.trajXSelect=uicontrol(gt.trajControlPanel,...          % Select the plot type
                'Style','popup','units','normalized',...
                'Position',[0.07,0.82,0.25,0.08],...'Position',[20,30,160,20],...
                'String',gt.trajXnames,...
                'Value',1,...
                'Callback',@gt.changeTrajectoryPlotVariable);
            
            uicontrol('Parent',gt.trajControlPanel,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','y',...
                'Units','normalized',...
                'Position',[.01 .68 .1 .1]);
            
            gt.trajYSelect=uicontrol(gt.trajControlPanel,...          % Select the plot type
                'Style','popup','units','normalized',...
                'Position',[0.07,0.72,0.25,0.08],...'Position',[20,30,160,20],...
                'String',gt.trajYnames,...
                'Value',1,...
                'Callback',@gt.changeTrajectoryPlotVariable);  
            
            
            uicontrol('Parent',gt.trajControlPanel,...
                'Style','text',...
                'HorizontalAlignment','left',...
                'String','z',...
                'Units','normalized',...
                'Position',[.01 .58 .1 .1]);
            
            gt.trajZSelect=uicontrol(gt.trajControlPanel,...          % Select the plot type
                'Style','popup','units','normalized',...
                'Position',[0.07,0.62,0.25,0.08],...'Position',[20,30,160,20],...
                'String',gt.trajZnames,...
                'Value',1,...
                'Callback',@gt.changeTrajectoryPlotVariable);  
            
            gt.trajTableDefault={[],[];[],[];[],[]}; %update the first time we plot... or populate with XPP values?
            
            gt.trajTable=uitable('Parent',gt.trajControlPanel,...
                'Units','normalized',...
                'Position',[0.32,0.5,0.68,0.5],...
                'RowName', {},...
                'ColumnName', {'min','max'},...
                'ColumnWidth',{'auto','auto'},...
                'ColumnFormat',{'numeric','numeric'},...
                'ColumnEditable', [true true],...%
                'Data',gt.trajTableDefault,...
                'CellEditCallback',@gt.editTrajPlotLims);
            
            %toggle buttons: snap to grid, subplots
            
            %parameter table (nclick, dt, t0, tf, nOut)            
            gt.numericsTableTraj=uitable('Parent',gt.trajControlPanel,...
                'Units','normalized',...
                'Position',[0.32,0.01,0.68,0.49],...
                'ColumnName', {},...
                'ColumnWidth','auto',...
                'ColumnFormat',{'numeric'},...
                'ColumnEditable', true,...
                'Data',gt.numericsTableDefaultTraj,...
                'RowName', gt.trajNumericsNames, ...
                'CellEditCallback',@gt.editTrajParams);
            
        end
        
    
%%%% Process key-presses from the different windows. Emulate XPP hotkeys where possible!
        
        
        %%%% keypress for grid window. 
        function processKeyPressGrid(gt,~,~)
            newKey=get(gcf,'currentkey');
            
            if gt.lastkeypress=='i'
                switch newKey
                    case 'g'
                        gt.integrate('go');
                    case 'l'
                        gt.integrate('last');
                    case 'r'
                        gt.integrate('random');
%                     case 's'
%                         gt.integrate('shift');
                    case 't'
                        gt.integrate('transient');
                end
                
            elseif gt.lastkeypress=='w'
                switch newKey
                    case 'f' %fit bounds to last solved values
                        gt.setBounds([gt.x(1),gt.y(1)],[gt.x(end),gt.y(end)]);
                    case 'o' %zoom out using rubber box
                        gt.zoomOut();
                    case 'z' %zoom in using rubber box
                        gt.zoomIn();
                end
            else
                switch newKey
%                     case 'c'
%                         gt.integrate('continue');
                    case 'add'
                        % + dz
                        gt.addToZ(gt.dz);
                        
                    case 'subtract'
                        % + dz
                        gt.addToZ(-gt.dz);
                    otherwise
                        gt.processCommonKeyPress(newKey);
                end
            end
                
            gt.lastkeypress=newKey;
        end
        
        %%%% keypress for control window. 
        function processKeyPressControl(gt,~,~)
            newKey=get(gcf,'currentkey');
            
            gt.processCommonKeyPress(newKey);
            
            gt.lastkeypress=newKey;
        end
        
        %%%% keypress for trajector windows. 
        function processKeyPressTrajectories(gt,~,~)
            newKey=get(gcf,'currentkey');
            
            if gt.lastkeypress=='i'
                switch newKey
                    case 'g'
                        gt.integrateTraj('go');
                    case 'l'
                        gt.integrateTraj('last');
                    case 'r'
                        gt.integrateTraj('random');
                    case 's'
                        gt.integrateTraj('shift');
                end
            elseif gt.lastkeypress=='w'
                switch newKey
                    case 'f'
                        %make axes tight for all trajectory axes present
                        gt.fitTrajAxis();
                        
%                     case 'o' 
                        %activate zoom out tool
%                     case 'z' 
                        %activate zoom in tool  %Doesn't work: for some
                        %reason dumps the 'z' to command window...
%                         zoom on 
%                         figure(gt.trajFig)
                end
            else
                switch newKey
                    case 'c'
                        gt.integrateTraj('continue');
                    case 'add'
                        % + dz
                        gt.addToZ(gt.dz);
                        
                    case 'subtract'
                        % + dz
                        gt.addToZ(-gt.dz);
                    otherwise
                        gt.processCommonKeyPress(newKey);
                end
            end
                
            gt.lastkeypress=newKey;
        end
        
        %%%% common keystrokes across all windows
        %is there a better choice for f1-f3?
        function processCommonKeyPress(gt, newKey)
            if gt.lastkeypress=='f'
                switch newKey
                    case 'q'
                        gt.my_closereq();
                end
            else
                switch newKey
                case 'f1'
                    if ~isvalid(gt.gridFig)
                        gt.buildGridFig();
                    end
                    figure(gt.gridFig)
                    
                case 'f2'
                    figure(gt.controlFig)
                    
                case 'f3'
                    if ~isvalid(gt.trajFig)
                        gt.buildTrajFig();
                    end
                    figure(gt.trajFig)
                    
%                     case 'p'
%                         uitable(gt.parTable);
%                     case 'r'
%                         uitable(gt.gridTable);
%                     case 'u'
%                         uitable(gt.numericsTable);
%                     case 'y'
%                         uitable(gt.y0Table);

                end
            end
        end
        
        function addToZ(gt,amountToAdd)
            gt.p(gt.zIx)=gt.p(gt.zIx)+amountToAdd;
            if abs(gt.p(gt.zIx))<1e-15
                gt.p(gt.zIx)=0;
            end
            gt.parTable.Data{gt.zIx,2}=gt.p(gt.zIx);
            gt.generateProblemData();
            set(gt.zValueText,'String',[gt.pNames{gt.zIx} '=' num2str(gt.p(gt.zIx))]);
        end
        
%%%% run integration of grids or trajectories
        
        
        %%%% grid integration handler function
        function integrate(gt, type)
            if ~(strcmp(type,'go')|strcmp(type,'random')|strcmp(type,'transient')) && ~gt.havePreviousGrid
                disp('No previous solution')
                return
            end
            
%             h = waitbar(0, 'Please wait...');
%             set(h, 'WindowStyle','modal', 'CloseRequestFcn','');

            tic
            switch type
%                     
%                 case 'continue'
%                     [tmpF, gt.Yf]=gt.clODEgrid.integrateContinue();
%                     gt.F=gt.fFunc(tmpF,gt.fNames);

                case 'go'
                    [tmpF, gt.Yf]=gt.clODEgrid.integrate();
                    [gt.F,~]=gt.fFunc(tmpF,gt.fNames);
                    gt.havePreviousGrid=true;
                    gt.wasTransient=false;
                    
                    
                    gt.f=gt.F(:,gt.featIx);
                    
                case 'last'
                    gt.Y0=gt.Yf;
                    gt.clODEgrid.setProblemData(gt.Y0(:),gt.Pars(:),gt.nPts);
                    [tmpF, gt.Yf]=gt.clODEgrid.integrate();
                    [gt.F,~]=gt.fFunc(tmpF,gt.fNames);
                    gt.wasTransient=false;
                    
                    gt.f=gt.F(:,gt.featIx);
                    
                    
                case 'random'
                    gt.randomizeY0();
                    [tmpF, gt.Yf]=gt.clODEgrid.integrate();
                    [gt.F,~]=gt.fFunc(tmpF,gt.fNames);
                    gt.havePreviousGrid=true;
                    gt.wasTransient=false;
                    
                    gt.f=gt.F(:,gt.featIx);
                    
%                 case 'shift'
%                     [tmpF, gt.Yf]=gt.clODEgrid.integrateShift();

                case 'transient'
                    if gt.havePreviousGrid
                        gt.Y0=gt.Yf;
                        gt.clODEgrid.setProblemData(gt.Y0(:),gt.Pars(:),gt.nPts);
                    end
                    gt.Yf=gt.clODEgrid.transient();
                    gt.havePreviousGrid=true;
                    gt.wasTransient=true;
                    
                    gt.f=gt.Yf(:,gt.obspars.varIx);

            end
            tElapsed=toc
            
%             waitbar(1, h,['Time Elapsed = ' num2str(tElapsed)]);
%             pause(.5)
%             delete(h)
            
            gt.ff=reshape(gt.f,[gt.ny,gt.nx]); 
            gt.plotGrid();
        end
        
        
        
        %%%% trajectory integration handler function
        function integrateTraj(gt, type)
            
            if ~gt.havePreviousTraj
                disp('No previous trajectory - click within the grid!')
                return
            end
            
%             h = waitbar(0, 'Please wait...');
%             set(h, 'WindowStyle','modal', 'CloseRequestFcn','');

            tic
            switch type
                case 'go'
                    
                    gt.TimesTraj=[gt.XPPopt.t0,gt.XPPopt.t0,diff(gt.TimesTraj([1,3]))];
                    gt.trajY0=gt.trajY0def;
                    
                    [t_,y_,yf_,dy_,aux_]=runTrajectories(gt);
                    gt.t=t_;
                    gt.traj=y_;
                    gt.trajYf=yf_;
                    gt.trajDY=dy_;
                    gt.trajAux=aux_;
                    
                case 'last'
                    
                    gt.TimesTraj=[gt.XPPopt.t0,gt.XPPopt.t0,diff(gt.TimesTraj([1,3]))];
                    gt.trajY0=gt.trajYf;
                    
                    [t_,y_,yf_,dy_,aux_]=runTrajectories(gt);
                    gt.t=t_;
                    gt.traj=y_;
                    gt.trajYf=yf_;
                    gt.trajDY=dy_;
                    gt.trajAux=aux_;
                    
                case 'random'
                    
                    gt.TimesTraj=[gt.XPPopt.t0,gt.XPPopt.t0,diff(gt.TimesTraj([1,3]))];
                    %[0,1], then rescale each column so that it lies in [y0lb,y0ub]
                    unitY0=rand(gt.nClick,gt.nVar); 
                    ranges=(gt.y0ub-gt.y0lb);
                    gt.trajY0=unitY0*diag( ranges(:)) + repmat(gt.y0lb(:)',gt.nClick,1);
                    
                    [t_,y_,yf_,dy_,aux_]=runTrajectories(gt);
                    gt.t=t_;
                    gt.traj=y_;
                    gt.trajYf=yf_;
                    gt.trajDY=dy_;
                    gt.trajAux=aux_;
                    
                case 'shift'
                    
                    gt.TimesTraj=[gt.trajTF,gt.trajTF,gt.trajTF+diff(gt.TimesTraj([1,3]))];
                    gt.trajY0=gt.trajYf;
                    
                    [t_,y_,yf_,dy_,aux_]=runTrajectories(gt);
                    gt.t=t_;
                    gt.traj=y_;
                    gt.trajYf=yf_;
                    gt.trajDY=dy_;
                    gt.trajAux=aux_;
                    
                case 'continue'
                    
                    gt.TimesTraj=[gt.trajTF,gt.trajTF,gt.trajTF+diff(gt.TimesTraj([1,3]))];
                    gt.trajY0=gt.trajYf;
                    
                    [t_,y_,yf_,dy_,aux_]=runTrajectories(gt);
                    gt.t=[gt.t,t_];
                    gt.traj=[gt.traj;y_];
                    gt.trajYf=yf_;
                    gt.trajDY=[gt.trajDY;dy_];
                    gt.trajAux=[gt.trajAux;aux_];
                    
%                     set(gt.trajSubPlotAx,'nextplot','add')
                    gt.trajAxisMode='auto';
%                     set(gt.trajSubPlotAx,'nextplot','replace')
            end
            
%             waitbar(1, h,['Time Elapsed = ' num2str(tElapsed)]);
%             pause(.5)
%             delete(h)
                    
            tElapsed=toc
            
            gt.plotTrajectory();

        end
        
        %%%% common steps for trajectories
        function [t,y,yf,dy,aux]=runTrajectories(gt)
            
            gt.clODEtraj.setSolverParams(gt.dtTraj, gt.TimesTraj);
            gt.clODEtraj.setProblemData(gt.trajY0(:),gt.trajPars(:),gt.nClick);
            [y,yf,dy,aux]=gt.clODEtraj.trajectory();
            t=linspace(gt.TimesTraj(1),gt.TimesTraj(end),length(y(:,1)));
            
        end
        
        
%%%% Functions that reset grid bounds + create problem data
        
        
        %%%% user specifies new grid bounds using mouse click and drag box
        function zoomIn(gt)
            axes(gt.gridAxis)
            set(gt.gridAxis,'buttondownfcn','')
            k=waitforbuttonpress;
            if k==0
                point1 = gt.gridAxis.CurrentPoint;    % location when mouse clicked
                rbbox;                 % rubberbox
                point2 = gt.gridAxis.CurrentPoint;    % location when mouse released
                point1 = point1(1,1:2);            % extract x and y
                point2 = point2(1,1:2);
                lb = min(point1,point2);           % calculate bounds
                ub = max(point1,point2);           % calculate bounds

                set(gt.gridAxis,'buttondownfcn',@gt.clickTrajectory)

                if all(ub-lb>0)
                    gt.setBounds(lb,ub);
                else
                    disp('Try selecting new bounds more slowly...')
                end
            else
                disp('Window Zoom-in canceled by keypress')
            end
        end
        
        
        %%%% user specifies new box into which old limits will be shrunk
        function zoomOut(gt)
            axes(gt.gridAxis)
            oldx=xlim();
            oldy=ylim();
            oldlb=[oldx(1),oldy(1)];
            oldub=[oldx(2),oldy(2)];
            oldrange=[diff(oldx),diff(oldy)];

            set(gt.gridAxis,'buttondownfcn','')
            k=waitforbuttonpress;
            if k==0
                point1 = gt.gridAxis.CurrentPoint;    % location when mouse clicked
                rbbox;                 % return figure units
                point2 = gt.gridAxis.CurrentPoint;    % location when mouse released. MUST BE WITHIN THE FIGURE!!!
                point1 = point1(1,1:2);            % extract x and y
                point2 = point2(1,1:2);
                selectlb = min(point1,point2);           % calculate bounds
                selectub = max(point1,point2);           % calculate bounds
                newrange = abs(point1-point2);       % and dimensions

                lambda=oldrange./newrange;
                lb=oldlb-lambda.*abs(oldlb-selectlb);
                ub=oldub+lambda.*abs(oldub-selectub);

                set(gt.gridAxis,'buttondownfcn',@gt.clickTrajectory)

                if all(newrange>0)
                    gt.setBounds(lb,ub);
                else
                    disp('Try selecting new bounds more slowly...')
                end
            else
                disp('Window Zoom-out canceled by keypress')
            end
        end
        
        
        %%%% helper function to set all bounds (for using zoom in, zoom out tools)
        function setBounds(gt,newlb,newub)

            %set new values in grid control panel
            gt.gridTable.Data{1,2}=newlb(1);
            gt.gridTable.Data{2,2}=newlb(2);
            gt.gridTable.Data{1,3}=newub(1);
            gt.gridTable.Data{2,3}=newub(2);

            %set new values in problem data and par/y0 control panel
            if gt.xispar
                gt.plb(gt.xIx)=newlb(1);
                gt.pub(gt.xIx)=newub(1);
                gt.parTable.Data{gt.xIx,3}=newlb(1);
                gt.parTable.Data{gt.xIx,4}=newub(1);
            else
                gt.y0lb(gt.xIx)=newlb(1);
                gt.y0ub(gt.xIx)=newub(1);
                gt.y0Table.Data{gt.xIx,3}=newlb(1);
                gt.y0Table.Data{gt.xIx,4}=newub(1);
            end

            if gt.yispar
                gt.plb(gt.yIx)=newlb(2);
                gt.pub(gt.yIx)=newub(2);
                gt.parTable.Data{gt.yIx,3}=newlb(2);
                gt.parTable.Data{gt.yIx,4}=newub(2);
            else
                gt.y0lb(gt.yIx)=newlb(2);
                gt.y0ub(gt.yIx)=newub(2);
                gt.y0Table.Data{gt.yIx,3}=newlb(2);
                gt.y0Table.Data{gt.yIx,4}=newub(2);
            end

            %shrink the axis to new bounds
            axis([newlb(1),newub(1),newlb(2),newub(2)])

            gt.generateProblemData();
            
            %remove old trajectory points if grid is changed - the
            %trajectories no longer represent valid points in grid window
            gt.havePreviousTraj=false;
            oldMarks=findobj('tag','clickMarks');
            if ~isempty(oldMarks)
                delete(oldMarks);
            end
        end
        

        %%%% update problem data and sync to clODE device
        function generateProblemData(gt)
            
            %update grid problem data
            gt.Pars=repmat(gt.p(:)',gt.nPts,1);
            gt.Y0=repmat(gt.y0(:)',gt.nPts,1);

            if gt.xispar
                gt.x=linspace(gt.plb(gt.xIx),gt.pub(gt.xIx),gt.nx);
            else
                gt.x=linspace(gt.y0lb(gt.xIx),gt.y0ub(gt.xIx),gt.nx);
            end
            if gt.yispar
                gt.y=linspace(gt.plb(gt.yIx),gt.pub(gt.yIx),gt.ny);
            else
                gt.y=linspace(gt.y0lb(gt.yIx),gt.y0ub(gt.yIx),gt.ny);
            end

            [gt.xx,gt.yy]=meshgrid(gt.x, gt.y);
            
            if gt.xispar
                gt.Pars(:,gt.xIx)=gt.xx(:);
            else
                gt.Y0(:,gt.xIx)=gt.xx(:);
            end
            if gt.yispar
                gt.Pars(:,gt.yIx)=gt.yy(:);
            else
                gt.Y0(:,gt.yIx)=gt.yy(:);
            end
            
            %sync to grid device
            gt.clODEgrid.setProblemData(gt.Y0(:),gt.Pars(:),gt.nPts);
            gt.gridIsDirty=true;
            
            %also update trajectory problem data (does not overwrite click
            %coordinates)
            if gt.havePreviousTraj
                if gt.xispar
                    clickx=gt.trajPars(:,gt.xIx);
                else
                    clickx=gt.trajY0(:,gt.xIx);
                end
                if gt.yispar
                    clicky=gt.trajPars(:,gt.yIx);
                else
                    clicky=gt.trajY0(:,gt.yIx);
                end

                gt.trajPars=repmat(gt.p(:)',gt.nClick,1);
                gt.trajY0=repmat(gt.y0(:)',gt.nClick,1);

                if gt.xispar
                    gt.trajPars(:,gt.xIx)=clickx(:);
                else
                    gt.trajY0(:,gt.xIx)=clickx(:);
                end
                if gt.yispar
                    gt.trajPars(:,gt.yIx)=clicky(:);
                else
                    gt.trajY0(:,gt.yIx)=clicky(:);
                end

                gt.clODEtraj.setProblemData(gt.trajY0(:),gt.trajPars(:),gt.nClick);
            end
        end
        
        
        
%%%% plotting functions        
        
        
        %%%% grid plotting function
        function plotGrid(gt)
            
            
            thisPlotType=gt.selectedPlotType;
            
            if length(unique(gt.f))==1
                thisPlotType='image'; 
            end
            
            
            if ~isvalid(gt.gridFig) || ~isvalid(gt.gridAxis)
                gt.buildGridFig(); %will ensure a valid gridAxis
            end
            
            axes(gt.gridAxis);
            delete(gt.gridPlotHandle)
            hold on
            
            switch thisPlotType
                case 'scatter'
                    
                    %parameters
                    marker='s'; %for aesthetic preference
                    markerarea=36; %markersize could be set set to a second feature, f2 or log(f2)
                    %markeredgecolor %none/black
                    %filled/not filled
                    
                    gt.gridPlotHandle=scatter(gt.xx(:),gt.yy(:),markerarea,gt.f,'filled',marker);
                    
%                 case 'surf'
%                     
%                     %parameters
%                     %edgecolor %none/black
%                     
%                     gt.gridPlotHandle=surf(gt.xx,gt.yy,gt.ff,'edgecolor','none');
%                     light('Position',[0 0 1],'Style','infinite');
%                     lighting phong
%                     axis vis3d
                    
                case 'image' 
                    
                    %no parameters
                    
                    gt.gridPlotHandle=imagesc(gt.x,gt.y,gt.ff);
                    set(gt.gridAxis,'ydir','normal')
                    
                case 'contourf'
                    
                    %parameters
                    leveltype='real'; %real/int
                    nLevels=10;
                    linespec='k-'; %-,--,:,-. + color
                    
                    switch leveltype
                        case 'real'
                            levels=nLevels;
                        case 'int'
                            levels=linspace(1,gt.cmax,gt.cmax-1)-0.0001;
                    end
                    
                    [C,gt.gridPlotHandle]=contourf(gt.xx,gt.yy,gt.ff,levels,linespec);
                    
                    %clabel %use instead of colorbar
            end
            
            if gt.wasTransient
                title([gt.vNames{gt.obspars.varIx} '_f'])
            else
                title(gt.fNames{gt.featIx})
%             title([gt.vNames{gt.clODEgrid.Obs.varIx} ': ' gt.fNames{gt.featIx}])
            end            

            if gt.xispar
%                 xlabel(gt.pNames{gt.xIx})
                xlabel(gt.gridAxis,gt.parTable.Data{gt.xIx,1})
            else
                xlabel(gt.v0Names{gt.xIx})
            end
            
            if gt.yispar
%                 ylabel(gt.pNames{gt.yIx})
                ylabel(gt.gridAxis,gt.parTable.Data{gt.yIx,1})
            else
                ylabel(gt.v0Names{gt.yIx})
            end
            
            colorbar
            
            axis tight
            
            if ~strcmpi(thisPlotType,'surf') %mouse click coordinates don't make sense in 3D plot
                set(gt.gridPlotHandle,'hittest','off') %allows clickthrough to axis for clickSolve
                set(gt.gridAxis,'buttondownfcn',@gt.clickTrajectory)
            end
            
            if strcmp(gt.cAxisMode,'auto')
                gt.fitCAxis();
            else
                set(gt.gridAxis,'clim',[gt.cmin,gt.cmax]);
            end
            
            oldMarks=findobj('tag','clickMarks');
            uistack(oldMarks,'top')
            gt.toggleShowParsMarker();
            gt.gridIsDirty=false;
        end
        
        
        %%%% mouse-click in main axis grabs (x,y) for trajectories
        function clickTrajectory(gt,~,~)
            
            gt.trajPars=repmat(gt.p(:)',gt.nClick,1);
            gt.trajY0=repmat(gt.y0(:)',gt.nClick,1);

            oldMarks=findobj('tag','clickMarks');
            if ~isempty(oldMarks)
                delete(oldMarks);
            end

                %if snap2grid toggled
                %else
                %   [px,py]=ginput(3);
                %end
                
            for i=1:gt.nClick
                [px,py]=ginput(1);
                gt.xClickVal(i)=px;
                gt.yClickVal(i)=py;
                
                dx=gt.x(2)-gt.x(1);
                dy=gt.y(2)-gt.y(1);
                r=find(abs(gt.x-px)<dx/2,1,'first');
                c=find(abs(gt.y-py)<dy/2,1,'first');
                
                if isempty(c) || isempty(r)
                    disp('Click must be inside the axis')
                    return
                end
                
                clickIx=sub2ind([gt.ny,gt.nx],c,r);
                
                %set current Y0 for trajectory
%                 gt.trajY0(i,:)=gt.Y0(sub2ind([gt.nx,gt.ny],c,r),:);
                gt.trajY0(i,:)=gt.Y0(clickIx,:);
                
                if gt.xispar
%                     gt.trajPars(i,gt.xIx)=gt.Pars(sub2ind([gt.nx,gt.ny],c,r),gt.xIx); %snap to grid  
                    gt.trajPars(i,gt.xIx)=gt.Pars(clickIx,gt.xIx); %snap to grid                     
                    gt.xClickVal(i)=gt.trajPars(i,gt.xIx);
                    
%                     gt.trajPars(i,gt.xIx)=px; %real value
                else
                    gt.trajY0(i,gt.xIx)=px;
                end

                if gt.yispar
%                     gt.trajPars(i,gt.yIx)=gt.Pars(sub2ind([gt.nx,gt.ny],c,r),gt.yIx); %snap to grid
                    gt.trajPars(i,gt.yIx)=gt.Pars(clickIx,gt.yIx); %snap to grid
                    gt.yClickVal(i)=gt.trajPars(i,gt.yIx);
                    
%                     gt.trajPars(i,gt.yIx)=py;
                else
                    gt.trajY0(i,gt.yIx)=py;
                end
                
                %print the features to command prompt
%                 gt.F(clickIx,:)
                
                line(gt.xClickVal(i),gt.yClickVal(i),'linewidth',1,...
                    'marker',gt.markers(i),'linestyle','none','color',gt.markerColor,...
                    'hittest','off','tag','clickMarks');
            end
            
            gt.trajY0def=gt.trajY0;
            gt.TimesTraj=[gt.XPPopt.t0,gt.XPPopt.t0,diff(gt.TimesTraj([1,3]))];
            
            gt.clODEtraj.setSolverParams(gt.dtTraj, gt.TimesTraj);
            gt.clODEtraj.setProblemData(gt.trajY0(:),gt.trajPars(:),gt.nClick);
            
            [gt.traj, gt.trajYf, gt.trajDY, gt.trajAux]=gt.clODEtraj.trajectory();
            gt.t=linspace(gt.TimesTraj(1),gt.TimesTraj(end),length(gt.traj(:,1)));
            gt.havePreviousTraj=true;
            
            gt.plotTrajectory();
            
        end
        
        
        %%%% plot trajectories
        function plotTrajectory(gt) 
            
            %%%% extract the relevant X and Y vectors for plotting
            
            %extract the presently selected plotting variables
            gt.trajTF=gt.t(end);
            ix=get(gt.trajXSelect,'value');
            iy=get(gt.trajYSelect,'value');
            iz=get(gt.trajZSelect,'value');
            
            
            trajPlotX=zeros(length(gt.t),gt.nClick);
            trajPlotY=zeros(length(gt.t),gt.nClick);
            trajPlotZ=zeros(length(gt.t),gt.nClick);
            for i=1:gt.nClick
                
                %X axis
                if ix==1 %TIME
                    trajPlotX(:,i)=gt.t;
%                     trajPlotX(:,i)=gt.t/1000;

                elseif ix<gt.nVar+2 %Variable
                    thisTraj=gt.traj(:,i:gt.nClick:end);
                    trajPlotX(:,i)=thisTraj(:,ix-1);
                    
                elseif ix<2*gt.nVar+2 %Derivative
                    thisTraj=gt.trajDY(:,i:gt.nClick:end);
                    trajPlotX(:,i)=thisTraj(:,ix-(gt.nVar+1));
                    
                else %Auxiliary quantity
                    thisTraj=gt.trajAux(:,i:gt.nClick:end);
                    trajPlotX(:,i)=thisTraj(:,ix-(2*gt.nVar+1));
                end

                %Y axis - should time be allowed here too?
                if iy<gt.nVar+1 %Variable
                    thisTraj=gt.traj(:,i:gt.nClick:end);
                    trajPlotY(:,i)=thisTraj(:,iy);
                    
                elseif iy<2*gt.nVar+1 %Derivative
                    thisTraj=gt.trajDY(:,i:gt.nClick:end);
                    trajPlotY(:,i)=thisTraj(:,iy-gt.nVar);
                    
                else %Auxiliary quantity
                    thisTraj=gt.trajAux(:,i:gt.nClick:end);
                    trajPlotY(:,i)=thisTraj(:,iy-(2*gt.nVar));
                end
                
                %Z axis - will plot in 2D if iz==1
                if iz==1
                elseif iz<gt.nVar+2 %Variable
                    thisTraj=gt.traj(:,i:gt.nClick:end);
                    trajPlotZ(:,i)=thisTraj(:,iz-1);
                    
                elseif iz<2*gt.nVar+2 %Derivative
                    thisTraj=gt.trajDY(:,i:gt.nClick:end);
                    trajPlotZ(:,i)=thisTraj(:,iz-(gt.nVar+1));
                    
                else %Auxiliary quantity
                    thisTraj=gt.trajAux(:,i:gt.nClick:end);
                    trajPlotZ(:,i)=thisTraj(:,iz-(2*gt.nVar+1));
                end
            end
            
            %%%% plot the data and format the plot
            
            if ~isvalid(gt.trajFig)
                gt.buildTrajFig();
            end
            figure(gt.trajFig)
            delete(findobj('tag','trajMarks'));
            
            %if subplotToggleOn use subplots
            %else use one axis with multiple lines
            
            Lims=[Inf(3,1),-Inf(3,1)];
            
            for i=1:gt.nClick
                
                gt.trajSubPlotAx(i)=subplot(gt.nClick,1,i);
                
                if iz==1
                    plot(trajPlotX(:,i),trajPlotY(:,i),'k')
                else
                    plot3(trajPlotX(:,i),trajPlotY(:,i),trajPlotZ(:,i),'k')
                    axis vis3d
                    grid on
                end
                
                %x-variable is TIME special items... ACTUALLY this is
                %plotting feature detector stuff. Maybe this can be toggled
                if ix==1 && iz==1
                    %indicator for the beginning of observation interval:
                    if gt.Times(2)<max(trajPlotX)
                        line(gt.Times(2)*[1,1],ylim(),'color',[0,0.5,0],'linestyle',':')
                    end

                    %if using plateau observer, plot the detector thresholds
                    if strcmpi(gt.observer,'plateau') && ix==1 && iy < gt.nVar+1
                        ymax=max(trajPlotY(:,i));
                        ymin=min(trajPlotY(:,i));
                        threshYup=ymin+gt.obspars.fractionYup*abs(ymax-ymin);
                        threshYdown=ymin+gt.obspars.fractionYdown*abs(ymax-ymin);
                        line([gt.TimesTraj(1),gt.TimesTraj(end)],[threshYup,threshYup],'color','r')
                        line([gt.TimesTraj(1),gt.TimesTraj(end)],[threshYdown,threshYdown],'color','b')
                    end
                end

                
                title([gt.gridXname '=' num2str(gt.xClickVal(i)) ', ' gt.gridYname '=' num2str(gt.yClickVal(i))]);
                
                if i~=gt.nClick && ix==1 && iz==1
                    set(gca,'xticklabel','')
                else
                    xlabel(gt.trajXnames{ix}); 
                end
                ylabel(gt.trajYnames{iy});
                
                if iz~=1
                    zlabel(gt.trajZnames{iz});
                end
                
                
                Xrange=[min(trajPlotX(:,i)),max(trajPlotX(:,i))];
                Yrange=[min(trajPlotY(:,i)),max(trajPlotY(:,i))];
                Zrange=[min(trajPlotZ(:,i)),max(trajPlotZ(:,i))];

                
                %set axis limits
                if strcmp(gt.trajAxisMode,'auto') %From table
                    XLIM=Xrange;
                    YLIM=Yrange;
                    ZLIM=Zrange;

                    if diff(XLIM)==0
                        if XLIM(1)==0
                            XLIM=[-1,1];
                        else
                            XLIM=XLIM+0.01*abs(XLIM(1))*[-1,1];
                        end
                    elseif any(isnan(XLIM))
                        XLIM=[0,1];
                    end
                    if diff(YLIM)==0
                        if YLIM(1)==0
                            YLIM=[-1,1];
                        else
                            YLIM=YLIM+0.01*abs(YLIM(1))*[-1,1];
                        end
                    elseif any(isnan(YLIM))
                        YLIM=[0,1];
                    end
                    if diff(ZLIM)==0
                        if ZLIM(1)==0
                            ZLIM=[-1,1];
                        else
                            ZLIM=ZLIM+0.01*abs(ZLIM(1))*[-1,1];
                        end
                    elseif any(isnan(ZLIM))
                        ZLIM=[0,1];
                    end
                else
                    XLIM=cell2mat(gt.trajTable.Data(1,:));
                    YLIM=cell2mat(gt.trajTable.Data(2,:));
                    ZLIM=cell2mat(gt.trajTable.Data(3,:));
                end
                
                set(gca,'xlim',XLIM,'ylim',YLIM);
                if iz~=1
                    set(gca,'zlim',ZLIM);
                end
                
                Lims(:,1)=min([XLIM(1);YLIM(1);ZLIM(1)],Lims(:,1)); 
                Lims(:,2)=max([XLIM(2);YLIM(2);ZLIM(2)],Lims(:,2));
                
                %add marker in upper right
                xMark=XLIM(2)-diff(XLIM)/50;
                yMark=YLIM(2)-diff(YLIM)/10;
                zMark=0;
                if iz~=1
                    zMark=ZLIM(2)-diff(ZLIM)/10;
                end
                
                line(xMark,yMark,zMark,'linewidth',1,...
                    'marker',gt.markers(i),'linestyle','none','color',gt.markerColor,...
                    'tag','trajMarks');
                
                box off
            end
            
            
            if strcmp(gt.trajAxisMode,'auto')
                gt.trajTable.Data(1,:)=num2cell(Lims(1,:));
                gt.trajTable.Data(2,:)=num2cell(Lims(2,:));
                gt.trajTable.Data(3,:)=num2cell(Lims(3,:));
            end
        end
        
        
        
        
        %UI callbacks
        
        %%% control window callbacks
        %TODO:  - error checking = lb<ub, ub>lb, nPts positive int
        %       - make sure no complex values allowed in table!
        %       - make a subroutine to update grid all data:
        %       gt.updateProblemData(newValue,argsToSpecifyWhatNeedsUpdating)???
        % ALT: use listeners??
        function editGridControls(gt, hTable, callbackdata)
            
            %row and column of edited cell
            row=callbackdata.Indices(1);
            col=callbackdata.Indices(2);
            NewData=callbackdata.NewData;
      
            % reject accidental entry of letters
            if isnan(NewData)
                hTable.Data{row,col}=callbackdata.PreviousData;
                return
            end
            
            % if accidentally press "i" while editing table, uitable
            % interprets as complex number. reject this
            if ~isreal(NewData)
                hTable.Data{row,col}=callbackdata.PreviousData;
                return
            end
            
            
            if col==1 %change identity of x or y
                if row==1
                    [gt.xispar,gt.xIx]=ismember(NewData,gt.pNames);
                    if gt.xispar
                        gt.gridTable.Data{row,2}=gt.plb(gt.xIx);
                        gt.gridTable.Data{row,3}=gt.pub(gt.xIx);
                        gt.gridXname=gt.pNames{gt.xIx};
                    else
                        [~,gt.xIx]=ismember(NewData,gt.v0Names);
                        gt.gridTable.Data{row,2}=gt.y0lb(gt.xIx);
                        gt.gridTable.Data{row,3}=gt.y0ub(gt.xIx);
                        gt.gridXname=gt.v0Names{gt.xIx};
                    end
                    
                else
                    [gt.yispar,gt.yIx]=ismember(NewData,gt.pNames);
                    if gt.yispar
                        gt.gridTable.Data{row,2}=gt.plb(gt.yIx);
                        gt.gridTable.Data{row,3}=gt.pub(gt.yIx);
                        gt.gridYname=gt.pNames{gt.yIx};
                    else
                        [~,gt.yIx]=ismember(NewData,gt.v0Names);
                        gt.gridTable.Data{row,2}=gt.y0lb(gt.yIx);
                        gt.gridTable.Data{row,3}=gt.y0ub(gt.yIx);
                        gt.gridYname=gt.v0Names{gt.yIx};
                    end
                end
%                 set(gt.gridAxis,'xlim',[gt.gridTable.Data{1,2},gt.gridTable.Data{1,3}],...
%                     'ylim',[gt.gridTable.Data{2,2},gt.gridTable.Data{2,3}]);
                 
                
            elseif col==2 %change in lower bound
                if row==1
                    if gt.xispar
                        gt.plb(gt.xIx)=NewData;
                        gt.parTable.Data{gt.xIx,3}=NewData;
                    else
                        gt.y0lb(gt.xIx)=NewData;
                        gt.y0Table.Data{gt.xIx,3}=NewData;
                    end
                else
                    if gt.yispar
                        gt.plb(gt.yIx)=NewData;
                        gt.parTable.Data{gt.yIx,3}=NewData;
                    else
                        gt.y0lb(gt.yIx)=NewData;
                        gt.y0Table.Data{gt.yIx,3}=NewData;
                    end
                end
                
%                 set(gt.gridAxis,'xlim',[gt.gridTable.Data{1,2},gt.gridTable.Data{1,3}],...
%                     'ylim',[gt.gridTable.Data{2,2},gt.gridTable.Data{2,3}]);
                
                if gt.xispar
                    xlabel(gt.gridAxis,gt.pNames{gt.xIx})
%                     xlabel(gt.gridAxis,gt.parTable.Data{gt.xIx,1})
                else
                    xlabel(gt.gridAxis,gt.v0Names{gt.xIx})
                end

                if gt.yispar
                    ylabel(gt.gridAxis,gt.pNames{gt.yIx})
                else
                    ylabel(gt.gridAxis,gt.v0Names{gt.yIx})
                end
            
                
            elseif col==3 %change in upper bound
                if row==1
                    if gt.xispar
                        gt.pub(gt.xIx)=NewData;
                        gt.parTable.Data{gt.xIx,4}=NewData;
                    else
                        gt.y0ub(gt.xIx)=NewData;
                        gt.y0Table.Data{gt.xIx,4}=NewData;
                    end
                else
                    if gt.yispar
                        gt.pub(gt.yIx)=NewData;
                        gt.parTable.Data{gt.yIx,4}=NewData;
                    else
                        gt.y0ub(gt.yIx)=NewData;
                        gt.y0Table.Data{gt.yIx,4}=NewData;
                    end
                end
%                 set(gt.gridAxis,'xlim',[gt.gridTable.Data{1,2},gt.gridTable.Data{1,3}],...
%                     'ylim',[gt.gridTable.Data{2,2},gt.gridTable.Data{2,3}]);
                
                
            elseif col==4 %change in number of points
                if row==1
                    gt.nx=NewData;
                    gt.nPts=gt.nx*gt.ny;
                else
                    gt.ny=NewData;
                    gt.nPts=gt.nx*gt.ny;
                end
                gt.havePreviousGrid=false;
            end
            
            set(gt.gridAxis,'xlim',[gt.gridTable.Data{1,2},gt.gridTable.Data{1,3}],...
                'ylim',[gt.gridTable.Data{2,2},gt.gridTable.Data{2,3}]);
                
            gt.generateProblemData();
            
            %remove old trajectory points if grid is changed - the
            %trajectories no longer represent valid points in grid window
            gt.havePreviousTraj=false;
            oldMarks=findobj('tag','clickMarks');
            if ~isempty(oldMarks)
                delete(oldMarks);
            end
            gt.toggleShowParsMarker();
        end
        
        
        function editParams(gt, hTable, callbackdata)
            
            %row and column of edited cell
            row=callbackdata.Indices(1);
            col=callbackdata.Indices(2);
            
            if isnan(callbackdata.NewData)
                hTable.Data{row,col}=callbackdata.PreviousData;
                return
            end
            
            if col==1
%                 gt.pNames{row}=callbackdata.NewData;
%                 tempix=find(strcmp(gt.gridTable.Data(:,1),callbackdata.NewData));
%                 if ~isempty(tempix)
%                     gt.gridTable.Data{tempix,1}=callbackdata.NewData;
%                 end
                
            elseif col==2
                gt.p(row)=callbackdata.NewData;
            elseif col==3
                gt.plb(row)=callbackdata.NewData;
                %update gridTable if x or y are parameters
                if gt.xispar && row==gt.xIx
                    gt.gridTable.Data{1,2}=callbackdata.NewData;
                end
                if gt.yispar && row==gt.yIx
                    gt.gridTable.Data{2,2}=callbackdata.NewData;
                end
            elseif col==4
                gt.pub(row)=callbackdata.NewData;
                if gt.xispar && row==gt.xIx
                    gt.gridTable.Data{1,3}=callbackdata.NewData;
                end
                if gt.yispar && row==gt.yIx
                    gt.gridTable.Data{2,3}=callbackdata.NewData;
                end
            end
            
            gt.generateProblemData();
            gt.toggleShowParsMarker();
        end
        
        
        function toggleShowParsMarker(gt,~,~)
            
            %get state of toggle
            toggleValue=get(gt.showParToggle,'value');
            
            %always remove old lines whether plotting new ones or not.
            oldLines=findobj('tag','pLine');
            if ~isempty(oldLines)
                delete(oldLines);
            end
            
            dx=gt.x(2)-gt.x(1);
            dy=gt.y(2)-gt.y(1);
            
            XLIM=xlim(gt.gridAxis);
            YLIM=ylim(gt.gridAxis);
            
            if toggleValue
                %vertical line (x-axis value)
                if gt.xispar
                    Xvert=gt.p(gt.xIx)*[1,1];
                else
                    Xvert=gt.y0(gt.xIx)*[1,1];
                end
                Yvert=YLIM;
                    
                %horizontal line (y-axis value)
                Xhoriz=XLIM;
                if gt.yispar
                    Yhoriz=gt.p(gt.yIx)*[1,1];
                else
                    Yhoriz=gt.y0(gt.yIx)*[1,1];
                end
                
                Xhoriz(1)=min(Xvert(1),Xhoriz(1))-dx;
                Xhoriz(2)=max(Xvert(2),Xhoriz(2))+dx;
                Yvert(1)=min(Yvert(1),Yhoriz(1))-dy;
                Yvert(2)=max(Yvert(2),Yhoriz(2))+dy;
                
                line(Xvert,Yvert,...
                    'parent',gt.gridAxis,'linewidth',1,...
                    'marker','none','linestyle','-','color','r',...
                    'hittest','off','tag','pLine');
                
                line(Xhoriz,Yhoriz,...
                    'parent',gt.gridAxis,'linewidth',1,...
                    'marker','none','linestyle','-','color','r',...
                    'hittest','off','tag','pLine');
                
                set(gt.gridAxis,'xlim',XLIM,'ylim',YLIM)
            end
        end
        
        function clickSetPars(gt,~,~)
            
            figure(gt.gridFig)
            [px,py]=ginput(1);
            
            %update gridTable if x or y are parameters
            if gt.xispar
                gt.parTable.Data{gt.xIx,2}=px;
                gt.p(gt.xIx)=px;
            else
                gt.y0Table.Data{gt.xIx,2}=px;
                gt.y0(gt.xIx)=px;
            end
            if gt.yispar
                gt.parTable.Data{gt.yIx,2}=py;
                gt.p(gt.yIx)=py;
            else
                gt.y0Table.Data{gt.yIx,2}=py;
                gt.y0(gt.yIx)=py;
            end
            gt.toggleShowParsMarker();
        end
        
        
        % Copy current parameter set into an XPP file (requires
        % ChangeXPPodeFile and packagePars4XPP)
        function savePars(gt,~,~)
            
            [fileToWrite,path]=uiputfile('*.ode','Save parameters to ODE file', gt.odeFilename);
            
            fileToWrite=fullfile(path,fileToWrite);
            
            if ~exist(fileToWrite,'file')
                copyfile(gt.odeFilename,fileToWrite)
            end
            
            ChangeXPPodeFile(fileToWrite,package4XPP(gt.pNames,gt.p));
        end
        
        
        function resetPars(gt,~,~)
            
            gt.p=[gt.parTableDefault{:,2}];
%             gt.plb=[gt.parTableDefault{:,3}];
%             gt.pub=[gt.parTableDefault{:,4}];
            
            tempParTableData=get(gt.parTable,'Data');
            tempParTableData(:,2)=gt.parTableDefault(:,2);
            set(gt.parTable,'Data',tempParTableData);
            
%             if gt.xispar
%                 gt.gridTable.Data{1,2}=gt.parTableDefault{gt.xIx,3};
%                 gt.gridTable.Data{1,3}=gt.parTableDefault{gt.xIx,4};
%             end
%             
%             if gt.yispar
%                 gt.gridTable.Data{2,2}=gt.parTableDefault{gt.yIx,3};
%                 gt.gridTable.Data{2,3}=gt.parTableDefault{gt.yIx,4};
%             end
            
%             set(gt.gridAxis,'xlim',[gt.gridTable.Data{1,2},gt.gridTable.Data{1,3}],...
%                 'ylim',[gt.gridTable.Data{2,2},gt.gridTable.Data{2,3}]);
            
%             if gt.xispar||gt.yispar
%                 %remove old trajectory points if grid is changed - the
%                 %trajectories no longer represent valid points in grid window
%                 gt.havePreviousTraj=false;
%                 oldMarks=findobj('tag','marks');
%                 if ~isempty(oldMarks)
%                     delete(oldMarks);
%                 end
%             end
            
            set(gt.zValueText,'String',[gt.pNames{gt.zIx} '=' num2str(gt.p(gt.zIx))]);
            gt.generateProblemData();
            gt.toggleShowParsMarker();
        end
        
        
        function editInitialData(gt, hTable, callbackdata)
            
            %row and column of edited cell
            row=callbackdata.Indices(1);
            col=callbackdata.Indices(2);
            
            if isnan(callbackdata.NewData)
                hTable.Data{row,col}=callbackdata.PreviousData;
                return
            end
            
            if col==1
                gt.v0Names{row}=callbackdata.NewData;
            elseif col==2
                gt.y0(row)=callbackdata.NewData;
            elseif col==3
                gt.y0lb(row)=callbackdata.NewData;
                %update gridTable if x or y are NOT parameters
                if ~gt.xispar && row==gt.xIx
                    gt.gridTable.Data{1,2}=callbackdata.NewData;
                end
                if ~gt.yispar && row==gt.yIx
                    gt.gridTable.Data{2,2}=callbackdata.NewData;
                end
            elseif col==4
                gt.y0ub(row)=callbackdata.NewData;
                if ~gt.xispar && row==gt.xIx
                    gt.gridTable.Data{1,3}=callbackdata.NewData;
                end
                if ~gt.yispar && row==gt.yIx
                    gt.gridTable.Data{2,3}=callbackdata.NewData;
                end
            end
            
            gt.generateProblemData();
        end
        
        
        function resetY0(gt,~,~)
            
            gt.y0=[gt.y0TableDefault{:,2}];
            gt.y0lb=[gt.y0TableDefault{:,3}];
            gt.y0ub=[gt.y0TableDefault{:,4}];
            
            set(gt.y0Table,'Data',gt.y0TableDefault);
            
            if ~gt.xispar
                gt.gridTable.Data{1,2}=gt.y0TableDefault{gt.xIx,3};
                gt.gridTable.Data{1,3}=gt.y0TableDefault{gt.xIx,4};
            end
            
            if ~gt.yispar
                gt.gridTable.Data{2,2}=gt.y0TableDefault{gt.yIx,3};
                gt.gridTable.Data{2,3}=gt.y0TableDefault{gt.yIx,4};
            end
            
            gt.generateProblemData();
        end
        
        
        function randomizeY0(gt,~,~)
            unitY0=rand(gt.nPts,gt.nVar); %[0,1]
            %rescale each column so that it lies in [y0lb,y0ub]
            ranges=(gt.y0ub-gt.y0lb);
            gt.Y0=unitY0*diag( ranges(:)) + repmat(gt.y0lb(:)',gt.nPts,1);
            
            %sync to grid device
            gt.clODEgrid.setProblemData(gt.Y0(:),gt.Pars(:),gt.nPts);
        end
        
        
        function editSolverParams(gt, hTable, callbackdata)
            row=callbackdata.Indices(1);
            
            if isnan(callbackdata.NewData)
                hTable.Data{row,1}=callbackdata.PreviousData;
                return
            end
            
            switch gt.gridNumericsNames{row}
                case 'dt'
                    gt.dt=callbackdata.NewData;
                case 't0'
                    gt.Times(1)=callbackdata.NewData;
                case 'trans'
                    gt.Times(2)=callbackdata.NewData;
                case 'total'
                    gt.Times(3)=callbackdata.NewData;
            end
            
            gt.clODEgrid.setSolverParams(gt.dt,gt.Times);
            gt.gridIsDirty=true;
        end
        
        
        function editTrajParams(gt, hTable, callbackdata)
            row=callbackdata.Indices(1);
            
            if isnan(callbackdata.NewData)
                hTable.Data{row,1}=callbackdata.PreviousData;
                return
            end
            
            switch gt.trajNumericsNames{row}
                case 'nClick'
                    gt.nClick=callbackdata.NewData;
                    gt.trajSubPlotAx=[];
                    gt.havePreviousTraj=false;
                case 'dt'
                    gt.dtTraj=callbackdata.NewData;
                case 't0'
                    gt.TimesTraj(1)=callbackdata.NewData;
                    gt.TimesTraj(2)=callbackdata.NewData;
                case 'total'
                    gt.TimesTraj(3)=gt.TimesTraj(1)+callbackdata.NewData;
                case 'nOut'
                    gt.nout=callbackdata.NewData;
            end
            
            gt.clODEtraj.setSolverParams(gt.dtTraj,gt.TimesTraj,gt.nout);
        end
        
        
        function editObserverParams(gt, hTable, callbackdata)
            
            row=callbackdata.Indices(1);
            col=callbackdata.Indices(2);
            
            if isnan(callbackdata.NewData)
                hTable.Data{row,col}=callbackdata.PreviousData;
                return
            end
            
            switch gt.obsparNames{row}
                case 'maxEventCount'
                    gt.obspars.maxEventCount=callbackdata.NewData;
                case 'minYamp'
                    gt.obspars.minYamp=callbackdata.NewData;
                case 'minDYamp'
                    gt.obspars.minDYamp=callbackdata.NewData;
                case 'minIMI'
                    gt.obspars.minIMI=callbackdata.NewData;
                case 'fractionYup'
                    gt.obspars.fractionYup=callbackdata.NewData;
                case 'fractionYdown'
                    gt.obspars.fractionYdown=callbackdata.NewData;
                case 'fractionDYup'
                    gt.obspars.fractionDYup=callbackdata.NewData;
                case 'fractionDYdown'
                    gt.obspars.fractionDYdown=callbackdata.NewData;
                case 'normTol'
                    gt.obspars.normTol=callbackdata.NewData;
            end
            gt.clODEgrid.setObserverParams(gt.obspars);
            gt.gridIsDirty=true;
            
        end
        
        
        function changeFeature(gt,callingObj,~)
            gt.featIx=get(callingObj,'value');
            
            if gt.havePreviousGrid
                gt.f=gt.F(:,gt.featIx);
                gt.ff=reshape(gt.f,[gt.nx,gt.ny]);
                gt.plotGrid();
                gt.fitCAxis(gt);
            end
        end
        
        
        function changeGridVariable(gt,callingObj,~)
            listNumber=get(callingObj,'value');
            gt.obspars.varIx=listNumber;
            
            gt.clODEgrid.setObserverParams(gt.obspars);
            gt.clODEtraj.setObserverParams(gt.obspars);
        end
        
        
        function changePlotType(gt,callingObj,~)
            listNumber=get(callingObj,'value');
            gt.selectedPlotType=gt.availablePlotTypes{listNumber};
            
            if gt.havePreviousGrid
                gt.plotGrid();
%                 gt.fitCAxis(gt);
            end
        end
        
        
        function editCmin(gt,callingObj,~)
            newval=str2double(get(callingObj,'string'));
            if newval<gt.cmax
                gt.cmin=newval;
                set(gt.gridAxis,'clim',[gt.cmin,gt.cmax]);
                gt.cAxisMode='manual';
            else
                set(gt.CAxisMinEdit,'string',num2str(gt.cmin));
            end
        end
        
        
        function editCmax(gt,callingObj,~)
            newval=str2double(get(callingObj,'string'));
            if newval>gt.cmin
                gt.cmax=newval;
                set(gt.gridAxis,'clim',[gt.cmin,gt.cmax]);
                gt.cAxisMode='manual';
            else
                set(gt.CAxisMaxEdit,'string',num2str(gt.cmax));
            end
        end
        
        
        function fitCAxis(gt,~,~)
            %expands the color limits to fit the data
            gt.cmin=min(gt.f);
            gt.cmax=max(gt.f); 
            
            if gt.cmax>gt.cmin
                set(gt.gridAxis,'clim',[gt.cmin,gt.cmax]);
            elseif gt.cmax==gt.cmin
                set(gt.gridAxis,'clim',[gt.cmin-1,gt.cmax+1]);
            end
            set(gt.CAxisMinEdit,'string',num2str(gt.cmin));
            set(gt.CAxisMaxEdit,'string',num2str(gt.cmax));
            gt.cAxisMode='auto';
        end
        
        
        function changeTrajectoryPlotVariable(gt,~,~)
            if gt.havePreviousTraj
                clf(gt.trajFig)
                gt.plotTrajectory();
            end
        end
        
        
        function fitTrajAxis(gt,~,~)
            gt.trajAxisMode='auto';
            if gt.havePreviousTraj
                clf(gt.trajFig)
                gt.plotTrajectory();
            end
        end
        
        
        function editTrajPlotLims(gt,~,~)
            gt.trajAxisMode='manual';
            XLIM=cell2mat(gt.trajTable.Data(1,:));
            YLIM=cell2mat(gt.trajTable.Data(2,:));
            set(gt.trajSubPlotAx,'xlim',XLIM,'ylim',YLIM);
        end
        
        
        function changeZpar(gt,~,~)
            gt.zIx=get(gt.zParSelect,'Value');
            if gt.p(gt.zIx)==0
                initialDZ=(gt.pub(gt.zIx)-gt.plb(gt.zIx))/10;
            else
                initialDZ=gt.p(gt.zIx)/10;
            end
            gt.dz=initialDZ;
            set(gt.dzEdit,'string',num2str(initialDZ))
            set(gt.zValueText,'String',[gt.pNames{gt.zIx} '=' num2str(gt.p(gt.zIx))]);
        end
        
        
        function editDZ(gt,~,callbackdata)
            oldDZ=gt.dz;
            NewData=str2double(get(gt.dzEdit,'String'));
%             zIx=get(gt.zParSelect,'Value');
            
            % reject accidental entry of letters
            if isnan(NewData)
                set(gt.dzEdit,'string',num2str(oldDZ))
                return
            end
            
            % if accidentally press "i" while editing table, uitable
            % interprets as complex number. reject this
            if ~isreal(NewData)
                set(gt.dzEdit,'string',num2str(oldDZ))
                return
            end
            
            if NewData>0 && isreal(NewData)
                gt.dz=NewData;
            end
        end
        
        function my_closereq(gt,~,~)
        % Close request function 
        % to display a question dialog box 
           selection = questdlg('Close Gridtool?',...
              'Close Request Function',...
              'Yes','No','Yes'); 
           switch selection,
              case 'Yes',
                 delete(gt.gridFig)
                 delete(gt.controlFig)
                 delete(gt.trajFig)
              case 'No'
              return 
           end
        end
        
    end
    
end