classdef gridtool < handle %matlab.mixin.SetGet
    %GRIDTOOL is a GUI for exploring the dependence of solutions of ODEs on
    %parameters, two-parameter planes at a time.  The ODEs are solved in
    %parallel using OpenCL via MEX interface to the c++ clODE library. The
    %interface is modeled after Bard Ermentrout's XPPAUT: ODE files that
    %run in XPP are parsed and converted to the required OpenCL code
    %automatically, and the interface shares some of the same keyboard
    %shortcuts.
    
    properties ( Dependent )
    end
    
    properties ( Access = private )
    end
    
    methods ( Access = public )
    end
    
    methods ( Access = private )
    end
    
end