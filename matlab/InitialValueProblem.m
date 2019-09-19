classdef InitialValueProblem

    
    properties ( Access = public )
    
        tscale=1
        tunit %for display
        
    end
    
        
    properties ( Access = private )
    end
    
    methods ( Access = public )
        
        function ivp=InitialValueProblem(varargin)
            %none - uigetfile to find ODE file
            %char - ODE filename
            %struct - parseODEfile output
        end
        
    end
    
    methods ( Access = private )
    end

end