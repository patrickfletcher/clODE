classdef clODEfeatures<clODE & matlab.mixin.SetGet
    % clODEfeatures(prob, stepper=rk4, observer=basic, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
   
    %TODO: this should have a method for default observer params
    
    %TODO: support post-processing functions for F - add amplitude etc
    %      [newF, newNames]=processFeatures(F, fNames);
    
    properties
        F
        observer='basic'
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
                observer='basic';
            end
%             observerInt=clODEfeatures.getObserverEnum(observer);
            
            if ~exist('mexFilename','var')
                mexFilename='clODEfeaturesmex';
            end
            obj@clODE(arg1, precision, selectedDevice, stepper, mexFilename, observer);
            
            obj.op=clODEfeatures.observerParams();
            obj.observerNames;
            obj.observer=observer;
        end
        
        %override initialize to include observerparams arg
        function initialize(obj, tspan, X0, P, sp, op)
            obj.cppmethod('initialize', tspan, X0(:), P(:), sp, op);
            obj.tspan=tspan;
            obj.X0=X0;
            obj.P=P;
            obj.sp=sp;
            obj.nPts=numel(X0)/obj.prob.nVar; 
            obj.op=op;
            obj.getNFeatures();
            obj.featureNames();
        end
        
        function setObserverPars(obj, op)
            if ~isequal(op,obj.op)
                obj.op=op;
                obj.cppmethod('setobserverpars', op);
            end
        end
        
        function set.observer(obj, newObserver)
            if ~strcmp(newObserver,obj.observer)
                if ismember(newObserver,obj.observerNames)
                    obj.observer=newObserver;
                    obj.cppmethod('setobserver', newObserver);
                    obj.featureNames();
                else
                    error(['undefined observer: ' newObserver]);
                end
            end
        end
        
        function initObserver(obj)
            obj.cppmethod('initobserver');
        end
        
        %overloads to fetch data if desired
        function features(obj, doInit)
            if ~exist('doInit','var')
                obj.cppmethod('features');
            else
                obj.cppmethod('features',doInit);
            end
        end
        
        
        function nf=getNFeatures(obj)
            obj.nFeatures=obj.cppmethod('getnfeatures');
            nf=obj.nFeatures;
        end
        
        function F=getF(obj, fix)
            F=obj.cppmethod('getf'); 
            F=reshape(F,obj.nPts,obj.nFeatures); %force column
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
    
     
    %static helper methods
    methods (Static=true)
        
        function op = observerParams()

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

