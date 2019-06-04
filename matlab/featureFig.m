classdef featureFig < handle
    %GRIDFIG class to display clODEfeatures on 2D domain
    %   simple class to encapsulate the plotting of features over 2D
    %   domain. main application: 2D grid
    
    properties (Access = private)
        clODEfeatsObj
        gridtoolObj
        
        FigH %figure handle
        AxesH %axis handle
        ImH %image handle from imagesc
        
        %listeners
        HLUpdateFigure
        HLLm
    end
    
    properties
        selectedFeature
    end
    
    methods
        function obj = featureFig(clODEfeatsObj, gridtoolObj, initialFeature)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.clODEfeatsObj = clODEfeatsObj;
            obj.gridtoolObj = gridtoolObj;
            obj.selectedFeature = initialFeature;
            obj.createLisn;
            
            obj.FigH=figure('NumberTitle','off',...
                'toolbar','figure',...
                'Name',['gridtool: ' clODEfeatsObj.ivp.name],...
                'Units','normalized',...
                'WindowKeyPressFcn',@gridtoolObj.processKeyPressGrid);
        end
        
        function createLisn(obj)
            obj.HLUpdateFigure = addlistener(obj.gridtoolObj,'UpdateGraph',...
                @(src,evnt)listenUpdateGraph(obj,src,evnt));
            obj.HLLm = addlistener(obj.gridtoolObj,'Lm','PostSet',...
                @(src,evnt)listenLm(obj,src,evnt));
        end
        
        function updateFigure(obj,inputArg)
        end
    end
end

