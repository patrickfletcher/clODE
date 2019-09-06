classdef gridfigure < Chart
    
    
    
    properties ( Dependent )
        % Chart x-data.
        XData
        % Chart y-data.
        YData
        % Color data for the scatter series.
        CData   
        % Axes x-label.
        XLabel
        % Axes y-label.
        YLabel
        % Axes title.
        Title
%         % Marker for clicked parameter points
%         MarkerSize
%         % Size data for the scatter series.
%         MarkerColor
    end
    
    properties ( Access = private )
        % Backing property for the x-data.
        XData_ = NaN;
        % Backing property for the y-data.
        YData_ = NaN;
        % Backing property for the y-data.
        CData_ = NaN;
        % image object
        Image
        % line object to show clicked points
        ClickedPoints
    end
    
    methods
    end
    
    methods ( Access = private )
    end
    
end