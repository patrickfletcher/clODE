classdef clODEfeatures < clODE
 
    %TODO: support post-processing functions for F - add amplitude etc
    %      [newF, newNames]=processFeatures(F, fNames);
    
    properties
        F
        observer
        op
        nFeatures
        oNames
        fNames
        
    %inherits from clODE:
%     prob
%     stepper
%     clSinglePrecision
%     cl_vendor
%     cl_deviceType
%         clBuilt=false
%         clInitialized=false
%     
%     nPts
%     P
%     X0
%     Xf
%     auxf
%     sp
%     tspan
    end
    
   
    methods
       
        function obj = clODEfeatures(arg1, precision, selectedDevice, stepper, observer, mexFilename)
            %hack to get correct mexfile for classes derived from this one.
            %I want the subclass to get all the methods contained here, but
            %it needs to use a mex function that unfortunately has to
            %repeat base class method dispatch code.
            
            if  ~exist('precision','var')||isempty(precision)
                precision=[]; %default handled in clODE.m
            end
            
            if  ~exist('selectedDevice','var')||isempty(selectedDevice)
                selectedDevice=[]; %default handled in clODE.m
            end
            
            if  ~exist('stepper','var')||isempty(stepper)
                stepper=[]; %default handled in clODE.m
            end
            
            if ~exist('observer','var')
                observer='basicall';
            end
            
            if ~exist('mexFilename','var')
                mexFilename='clODEfeaturesmex';
            end
            obj@clODE(arg1, precision, selectedDevice, stepper, mexFilename, observer);
            
            obj.observerNames();
            obj.observer=observer;
            obj.featureNames();
            obj.getNFeatures();
            obj.op=clODEfeatures.defaultObserverParams();%default observer params (device transfer during init)
        end
        
        %override initialize to include observerparams arg
        function initialize(obj, tspan, X0, P, sp, op)
            if ~exist('tspan','var') %no input args: use stored values
                tspan=obj.tspan; X0=obj.X0; P=obj.P; sp=obj.sp; op=obj.op;
            end
            obj.cppmethod('initialize', tspan, X0(:), P(:), sp, op);
            obj.tspan=tspan;
            obj.X0=X0;
            obj.P=P;
            obj.sp=sp;
            obj.nPts=numel(X0)/obj.prob.nVar; 
            obj.op=op;
            obj.featureNames();
            obj.getNFeatures();
            obj.clInitialized=true;
        end
        
        function setObserverPars(obj, op)
            if ~exist('op','var') %no input args: use stored values
                op=obj.op;
            end
            obj.op=op;
            obj.cppmethod('setobserverpars', op);
        end
        
        function setObserver(obj, newObserver)
            if ismember(newObserver,obj.observerNames)
                obj.observer=newObserver;
                obj.cppmethod('setobserver', newObserver);
                obj.featureNames();
                obj.getNFeatures();
                obj.clBuilt=false;
            else
                error(['undefined observer: ' newObserver]);
            end
        end
        
        function initObserver(obj)
            obj.cppmethod('initobserver');
        end
        
        function F=features(obj, doInit)
            if ~exist('doInit','var')
                obj.cppmethod('features');
            else
                obj.cppmethod('features',doInit);
            end
            if nargout==1 %overloads to fetch data if desired
                F=obj.getF();
            end
        end
        
        
        function nf=getNFeatures(obj)
            obj.nFeatures=obj.cppmethod('getnfeatures');
            nf=obj.nFeatures;
        end
        
        function F=getF(obj, fix)
            F=obj.cppmethod('getf'); 
            F=reshape(F,obj.nPts,obj.nFeatures);
            obj.F=F;
            if nargin==2 
                %return argument is just the subset fix of features 
                F=F(:,fix);
            end
        end
        
        function fNames=featureNames(obj)
            fNames=obj.cppmethod('getfeaturenames');
            obj.fNames=fNames;
        end
        
        function oNames=observerNames(obj)
            oNames=obj.cppmethod('getobservernames');
            obj.oNames=oNames;
        end
        
    end
    
    %ui public methods
    methods (Access=public)
        function hop=uisetOP(obj,parent)
            if ~exist('parent','var')||isempty(parent)
                parent=uifigure();
            end
            %TODO: convert observer pars to name/value struct..?
            sptable=table;
            sptable.name=fieldnames(obj.op);
            sptable.value=cell2mat(struct2cell(obj.op));
            hop=uitable(parent,'Data',sptable);
            hop.ColumnName = sptable.Properties.VariableNames;
%             hop.ColumnName = {'name'; 'value'};
            hop.ColumnWidth = {'auto', 'fit'};
            hop.RowName = {};
            hop.ColumnSortable = [false false];
            hop.ColumnEditable = [false true];
%             hop.Position = position;
            hop.CellEditCallback=@updateOP;
            
            function updateOP(src,event)
                thisfield=hop.DisplayData.name{event.Indices(1)};
                if isnumeric(event.NewData) && event.NewData>0
                    obj.op.(thisfield)=event.NewData;
                end
            end
        end
    end
     
    %static helper methods
    methods (Static=true)
        
        function op = defaultObserverParams()
            op.eVarIx=1; %not implemented
            op.fVarIx=1; %feature detection variable
            op.maxEventCount=5000; %stops if this many events found
            op.minXamp=0;  %don't record oscillation features if units of variable fVarIx
            op.minIMI=0; %not implemented
            op.nHoodRadius=0.25;  %size of neighborhood
            op.xUpThresh=0.3;  %not implemented
            op.xDownThresh=0.2; %selecting neighborhood centerpoint: first time fVarIx drops below this fraction of its amplitude {nhood2}
            op.dxUpThresh=0;  %not implemented
            op.dxDownThresh=0; %not implemented
            op.eps_dx=1e-6; %for checking for min/max
        end
        
        
    end
    
end

