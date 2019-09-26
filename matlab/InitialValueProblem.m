classdef InitialValueProblem < matlab.mixin.SetGet
    %A simple class to represent an initial value problem, parsed from an
    %XPP ode file.
    %
    % - provides a simple GUI to modify parameters and initial conditions
    % - 
    
    properties
    
        rhsfun %filname or matlab function that computes the vector field
        
        p
        x0
        
        tscale=1
        tunit='' %for display
        
    end
    
    properties ( Access = private )
        
    end
    
    methods
        
        function ivp=InitialValueProblem(varargin)
            %none - uigetfile to find ODE file
            %char - ODE filename
            %struct - parseODEfile output
        end
        
    end
    
    methods ( Access = private )
    end

end