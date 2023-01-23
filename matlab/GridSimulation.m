% a light wrapper class to represent a clODE 2-parameter grid simulation
% - basically want to provide grid construction/simulation convenience
% functions, expose clODEfeatures properties/methods with simple api.

%grid parameters can be par or init
%features can be var or aux

% TODO: A lot of the functionality here should really be built into the
% base clODE classes
% - the main C++ interface should hide conversion from table to problem
% struct expected by C++
classdef GridSimulation < handle

    properties (Access = public)
        name %a short name, defaulting to the odefile without extension

        tscale=1
        tlabel='Time'

        nX = 32
        nY = 32

        %table with all the allowed grid variables (her for now for easy setting..)
        gridvars

        grid=table('Size',[2,7], 'VariableTypes',...
            {'string','double','double','double','string','categorical','double'},...
            'VariableNames',{'name','value','lb','ub','label','type','ix'},...
            'RowNames',{'x','y'}) %table specifying grid x, y info
    end

    properties (SetAccess = private, GetAccess=public)
        %%%%%%%%%%%
        % set only during construction
        odefile %xpp ODE file - change would trigger entire refresh
        devices
        grid_device %must be in "devices" - trigger setOpenCL
        clo_g %clODE grid feature solver
        precision = 'single'
        stepper = 'dopri5'
        observer = 'basicall'

        % a dictionary listing features/computed feature names/functions
        featureMap

        GX %grid x, y as NDGRIDs
        GY

        XF %final state from previous sim

        %%%%%%%%%%%
        % change only by setter functions: want this to propagate changes
        tspan
        solverPars struct = clODEfeatures.defaultSolverParams()
        observerPars struct = clODEfeatures.defaultObserverParams()
    end

    % these are computed on demand from other properties
    properties (Dependent = true)
        nPts %

        %1D grids for x, y coords - computed from xvar, yvar
        gridx
        gridy

        % default parameter and initial condition points
        p0 % default parameter point
        x0 % default ic point

        % initial parameter and initial condition grids for simulation
        P0 %currently set parameter, IC matrices for clODE
        X0

        featureNames
    end

    methods
        function self = GridSimulation(odefile, opts)
            arguments
                odefile=[]
                opts.stepper='dopri5'
                opts.observer='basicall'
                opts.precision='single'
                opts.device_type="GPU"
                opts.nX=32
                opts.nY=[]
                opts.tscale=1
            end
            self.devices=queryOpenCL();
            self.grid_device=find({self.devices(:).type}==opts.device_type,1,'first');
            self.clo_g=clODEfeatures(odefile, opts.precision, self.grid_device, opts.stepper, opts.observer);

            if self.clo_g.prob.nPar+self.clo_g.prob.nVar<2
                error('Not enough parameters and variables for a 2D grid!')
            end

            %tspan
            t0=self.clo_g.prob.opt.t0;
            tf=self.clo_g.prob.opt.total;
            self.tspan=[t0,tf];
            self.clo_g.settspan([t0,tf]);
            self.nX=opts.nX;
            if isempty(opts.nY), self.nY=self.nX; end
            self.tscale=opts.tscale;

            %valid grid variables: parameters and/or initial conditions
            pars = self.clo_g.prob.par_tab(:,["name","value","lb","ub"]);
            pars.label = pars.name;
            pars.type = repmat("par",size(pars,1),1);
            pars.ix = (1:size(pars,1))';

            vars = self.clo_g.prob.var_tab(:,["name","value","lb","ub"]);
            vars.name = vars.name+"0";
            vars.label = vars.name;
            vars.type = repmat("ic",size(vars,1),1);
            vars.ix = (1:size(vars,1))';

            newGridvars=[pars;vars];
            newGridvars.type=categorical(newGridvars.type);
            newGridvars.Properties.RowNames=newGridvars.name;

            const_bounds=newGridvars.lb==newGridvars.ub;
            const_vals=newGridvars.value(const_bounds);
            newGridvars.lb(const_bounds)=const_vals-abs(const_vals)*0.1;
            newGridvars.ub(const_bounds)=const_vals+abs(const_vals)*0.1;

            self.gridvars = newGridvars;

            self.clo_g.buildCL();
            self.getFeatureMap();

            % user must then:
            % - set the grid names
            % - optionally change sp, op
            % - initialize the P0 and X0

            %             %set grid to first two names ???
            %             self.grid.name=newGridvars.name(1:2);
            %             self.gridvars=newGridvars; %triggers postSet listener
            %
            %             self.clo_g.Xf=self.X0;
            %             self.makeGridData();
            %             self.clo_g.initialize();
        end

        % change grid variable - optionally also set extra info
        function setGridVar(self, var, name, label, lb, ub, N)
            arguments
                self
                var
                name
                label = []
                lb = []
                ub = []
                N = []
            end

            for i=1:length(var)
                newvarinfo = self.gridvars(name(i),:);
                newvarinfo.Properties.RowNames = var(i);
                self.grid(var(i),:) = newvarinfo;
            end
        end

        function initialize(self, opts)
            arguments
                self
                opts.gridVars = []
                opts.ictype = "point"
                opts.newsp = []
                opts.newop = []
            end
            if isempty(opts.gridVars) && ismissing(self.grid("x","name")) && ismissing(self.grid("y","name"))
                self.setGridVar("x",self.gridvars.name(1));
                self.setGridVar("y",self.gridvars.name(2));
            end
            self.makeProblemData()
            self.clo_g.initialize()
        end

        function setTspan(self, newTspan)
            self.tspan=newTspan;
            self.clo_g.settspan(newTspan);
        end

        function makeProblemData(self, ictype)
            arguments
                self
                ictype = "point"
            end
            newP = repmat(self.p0,self.nPts,1);
            [self.GX,self.GY]=ndgrid(self.gridx, self.gridy);

            if self.grid{'x','type'}=="par"
                newP(:,self.grid{'x','ix'})=self.GX(:);
            end
            if self.grid{'y','type'}=="par"
                newP(:,self.grid{'y','ix'})=self.GY(:);
            end

            newX0 = self.makeInitialCondition(ictype);

            self.clo_g.setProblemData(newX0, newP);
        end

        function newX0=makeInitialCondition(self, type)
            %assumes makeGridData was run
            switch type
                case "last"
                    newX0=self.XF;
                case "point"
                    newX0=repmat(self.x0,self.nPts,1);
                case "random"
                    x0lb=self.gridvars.lb(self.gridvars.type=="ic")';
                    x0ub=self.gridvars.ub(self.gridvars.type=="ic")';
                    newX0=x0lb+rand(self.clo_g.nPts,length(x0lb)).*(x0ub-x0lb);
            end
            %make sure to overwrite any IC variable that is in the grid
            if self.grid{'x','type'}=="ic"
                newX0(:,self.grid{'x','ix'})=self.GX(:);
            end
            if self.grid{'y','type'}=="ic"
                newX0(:,self.grid{'y','ix'})=self.GY(:);
            end
        end


        % get Dependent values
        function nPts=get.nPts(self)
            nPts=self.nX * self.nY;
        end
        function p0=get.p0(self)
            p0=self.gridvars.value(self.gridvars.type=="par")';
        end
        function x0=get.x0(self)
            x0=self.gridvars.value(self.gridvars.type=="ic")';
        end
        function P0=get.P0(self)
            P0=self.clo_g.P;
        end
        function X0=get.X0(self)
            X0=self.clo_g.X0;
        end
        function XF=get.XF(self)
            XF=self.clo_g.XF;
        end
        function gridx=get.gridx(self)
            gridx=linspace(self.grid{"x","lb"},self.grid{"x","ub"},self.nX);
        end
        function gridy=get.gridy(self)
            gridy=linspace(self.grid{"y","lb"},self.grid{"y","ub"},self.nY);
        end
        function featureNames=get.featureNames(self)
            featureNames=string(self.featureMap.keys());
        end

        function getFeatureMap(self)
            self.clo_g.getNFeatures(); %update features (nVar may affect)
            fNames=self.clo_g.featureNames();
            f=containers.Map; %dict {name:func}

            %computed values
            ampvars=fNames(startsWith(fNames,'max')); %same as trajvars
            ampvars=split(ampvars);
            ampvars=ampvars(:,2);
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
            self.featureMap=f;
        end


        function integrate(self, action)
            %make sure all device variables are populated
            if ~self.clo_g.clInitialized
                self.initialize()
            end

            tic
            switch action
                case 'continue'
                    self.clo_g.shiftX0(); %device X0<-XF
                    self.clo_g.features(0);
                    self.clo_g.getF();

                case 'go'
                    self.clo_g.shiftX0();
                    self.clo_g.features(1);
                    self.clo_g.getF();

                    %                 case 'shift'
                case 'point'
                    newX0 = self.makeInitialCondition("point");
                    self.clo_g.setX0(newX0);
                    self.clo_g.transient();

                case 'pointgo'
                    newX0 = self.makeInitialCondition("point");
                    self.clo_g.setX0(newX0);
                    self.clo_g.transient();
                    self.clo_g.shiftX0();
                    self.clo_g.features(1);
                    self.clo_g.getF();

                case 'random'
                    newX0 = self.makeInitialCondition("random");
                    self.clo_g.setX0(newX0);
                    self.clo_g.transient();

                case 'transient'
                    self.clo_g.shiftX0();
                    self.clo_g.transient();
            end
            self.XF=self.clo_g.getXf();
            toc
        end

        function F = getFeature(self, featureName, fScale)
            arguments
                self
                featureName
                fScale = 1
            end
            func = self.featureMap(featureName);
            F = func(self.clo_g.F)/fScale;
            F = reshape(F,[self.nX,self.nY])';
        end

        function hIm = plotFeature(self, featureName, opts)
            arguments
                self
                featureName
                opts.featureLabel = []
                opts.fScale = 1
                opts.fig = []
                opts.ax = []
                opts.colormap = [0.9*[1,1,1];turbo]
                opts.cbgap = 0.025
                opts.cbdims = [0.3,0.05]
                opts.cbloc = "east"
                opts.cbjust = "midhi"
                opts.FontSize = 8
                opts.transform = @(F)F
            end

            F = self.getFeature(featureName, opts.fScale);
            F = opts.transform(F);

            F(F>1e10|F<-1e10) = nan;

%             fh = opts.fig;
%             if isempty(fh)
%                 fh=figure();
%             end
            ax = opts.ax;
            if isempty(ax)
                ax = gca;
            end
            featureLabel=opts.featureLabel;
            if isempty(opts.featureLabel)
                featureLabel=featureName;
            end

            hi=imagesc(ax, self.gridx, self.gridy, F);
            hi.HitTest='off';
            ax.YDir='normal';
            xlabel(self.grid.label(1));
            ylabel(self.grid.label(2));
%             axis(ax,'square')
            
%             hcb=colorbar('northoutside');
            [hcb, cblab] = makeCB(ax, ax.Position,opts.cbgap,opts.cbdims,opts.cbloc,opts.cbjust,featureLabel);
%             title(hcb,featureLabel,'Interpreter','none')

            colormap(ax,opts.colormap)

            ax.FontSize = opts.FontSize;
            hcb.FontSize = opts.FontSize;

            hIm.ax = ax;
            hIm.hi = hi;
            hIm.hcb = hcb;
            hIm.cblab = cblab;
        end

    end
end