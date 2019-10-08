classdef clODEfeatures<clODE
    % clODEfeatures(prob, stepper=rk4, observer=basic, clSinglePrecision=true, cl_vendor=any, cl_deviceType=default)
   
    %TODO: this should have a method for default observer params
    
    %TODO: support post-processing functions for F - add amplitude etc
    %      [newF, newNames]=processFeatures(F, fNames);
    
    properties
        fNames
        F
        observer
        op
        nFeatures
        
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
                observer='localmax';
            end
            observerInt=clODEfeatures.getObserverEnum(observer);
            
            if ~exist('mexFilename','var')
                mexFilename='clODEfeaturesmex';
            end
            obj@clODE(arg1, precision, selectedDevice, stepper, mexFilename, observerInt);
            
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
            obj.getFeatureNames();
        end
        
        function setObserverPars(obj, op)
            obj.op=op;
            obj.cppmethod('setobserverpars', op);
        end
        
        function set.observer(obj, newObserver)
            observerint=clODEfeatures.getObserverEnum(newObserver);
            obj.observer=newObserver;
            obj.cppmethod('setobserver', observerint);
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
        
        function F=getF(obj)
            F=obj.cppmethod('getf'); 
            F=reshape(F,obj.nPts,obj.nFeatures); %force column
            obj.F=F;
        end
        
        
        
        function fNames=getFeatureNames(obj)
            fNames={};
            switch lower(obj.observer)
                case {'basic'}
                    fNames{end+1}=['max ' obj.prob.varNames{obj.op.fVarIx}];
                    fNames{end+1}=['min ' obj.prob.varNames{obj.op.fVarIx}];
                    fNames{end+1}=['mean ' obj.prob.varNames{obj.op.fVarIx}];
                    fNames{end+1}=['max d' obj.prob.varNames{obj.op.fVarIx} '/dt'];
                    fNames{end+1}=['min d' obj.prob.varNames{obj.op.fVarIx} '/dt'];
                    fNames{end+1}='step count';
                    
                case {'basicallvar','basicall'}
                    for i=1:obj.prob.nVar
                        fNames{end+1}=['max ' obj.prob.varNames{i}];
                        fNames{end+1}=['min ' obj.prob.varNames{i}];
                        fNames{end+1}=['mean ' obj.prob.varNames{i}];
                        fNames{end+1}=['max d' obj.prob.varNames{i} '/dt'];
                        fNames{end+1}=['min d' obj.prob.varNames{i} '/dt'];
                    end
                    
                    for i=1:obj.prob.nAux
                        fNames{end+1}=['max ' obj.prob.auxNames{i}];
                        fNames{end+1}=['min ' obj.prob.auxNames{i}];
                        fNames{end+1}=['mean ' obj.prob.auxNames{i}];
                    end
                    fNames{end+1}='step count';
                    
                case {'localmax'}
                    fNames={...
                        'amplitude',...
                        'max IMI',...
                        'min IMI',...
                        'mean IMI',...
                        'max local amp',...
                        'min local amp',...
                        'mean local amp',...
                        ['max ' obj.prob.varNames{obj.op.fVarIx} '_{max}'],...
                        ['min ' obj.prob.varNames{obj.op.fVarIx} '_{max}'],...
                        ['mean ' obj.prob.varNames{obj.op.fVarIx} '_{max}'],...
                        ['max ' obj.prob.varNames{obj.op.fVarIx} '_{min}'],...
                        ['min ' obj.prob.varNames{obj.op.fVarIx} '_{min}'],...
                        ['mean ' obj.prob.varNames{obj.op.fVarIx} '_{min}'],...
                        'dx max',...
                        'dx min',...
                        ['mean ' obj.prob.varNames{obj.op.fVarIx}],...
                        'event count',...
                        'step count',...
                        };
%                         'max tMaxMin',...
%                         'min tMaxMin',...
%                         'mean tMaxMin',...

                case {'section1'}
                    
                case {'neighborhood1','nhood1'}
                    
                case {'section2'}
                    
                case {'neighborhood2','nhood2'}
                    fNames={...
                        'max period',...
                        'min period',...
                        'mean period',...
                        'max nMaxima',...
                        'min nMaxima',...
                        'mean nMaxima',...
                        'dx max',...
                        'dx min',...
                        ['mean ' obj.prob.varNames{obj.op.fVarIx}],...
                        'event count',...
                        'step count',...
                        };
                    
                otherwise %'basic'
                    fNames{end+1}=['max ' obj.prob.varNames{obj.op.fVarIx}];
                    fNames{end+1}=['min ' obj.prob.varNames{obj.op.fVarIx}];
                    fNames{end+1}=['mean ' obj.prob.varNames{obj.op.fVarIx}];
                    fNames{end+1}=['max d' obj.prob.varNames{obj.op.fVarIx} '/dt'];
                    fNames{end+1}=['min d' obj.prob.varNames{obj.op.fVarIx} '/dt'];
                    fNames{end+1}='nSteps';
                    
            end
            obj.fNames=fNames;
        end
        
    end
    
    
    
    %static helper methods
    methods (Static=true)
    
        function observerint=getObserverEnum(observername)
            switch lower(observername)
                case {'basic'}
                    observerint=0;
                case {'basicallvar','basicall'}
                    observerint=1;
                case {'localmax'}
                    observerint=2;
                case {'section1'}
                    observerint=3;
                case {'neighborhood1','nhood1'}
                    observerint=4;
                case {'section2'}
                    observerint=5;
                case {'neighborhood2','nhood2'}
                    observerint=6;
                otherwise
                    warning('Unrecognized observer method name. Using default: basic')
                    observerint=0;
            end
        end
        
    end
    
end

