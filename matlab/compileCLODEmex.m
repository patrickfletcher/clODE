% Script to compile clODE mex interface, so it can be used from Matlab.
% Make sure to set the correct paths to files in the Configuration block!
close all
clear

%TODO: auto detect paths for some common setups? e.g. CUDA

thisdir = pwd;

% %%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ismac % Code to run on Mac plaform
    opencl_include_dir = pwd; %cl.hpp for OpenCL C++ bindings
    opencl_lib_dir = ''; %leave empty; taken care of by the -framework option
    libopencl='';
    compflags='COMPFLAGS="$COMPFLAGS -std=c++11 -framework OpenCL"';
    ldflags='LDFLAGS="$LDFLAGS -framework OpenCL"';
    
elseif isunix % Code to run on Linux plaform
%     opencl_include_dir = pwd; %cl.hpp
    opencl_include_dir = '/usr/include'; %cl.hpp
    opencl_lib_dir = '/usr/local/lib'; %libOpenCL.so
    libopencl='-lOpenCL';
    compflags='COMPFLAGS="$COMPFLAGS -std=c++11"';
%     compflags='COMPFLAGS="$COMPFLAGS -W -Wall -Werror -ansi -pedantic"';
    ldflags='';
    opencl_lib_dir=['-L',opencl_lib_dir];
    
elseif ispc % Code to run on Windows platform 
%     opencl_include_dir = pwd; %cl.hpp
    opencl_include_dir = 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.1/include'; %cl.hpp
    opencl_lib_dir = 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.1/lib/x64';  %OpenCL.lib
    libopencl='-lOpenCL';
    compflags='COMPFLAGS="$COMPFLAGS /std:c++latest"';
%     compflags='COMPFLAGS="$COMPFLAGS -Wall"';
    ldflags='';
    opencl_lib_dir=['-L',opencl_lib_dir];
    
else
    disp('Cannot recognize platform')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd ../clode/cpp
clode_path=[pwd filesep]; 
clode_path=strrep(clode_path,'\','/')

spdlog_path = 'D:/GitHub/spdlog/include';

cd(thisdir)

debugchar='';
% debugchar='-g';

verbosechar='';
% verbosechar='-v';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compile commands:

%% queryOpenCL

mex('queryOpenCL.cpp',[clode_path,'OpenCLResource.cpp'],...
    debugchar,verbosechar,compflags,...
    ['-I' spdlog_path], ['-I' clode_path], ['-I' opencl_include_dir],...
    ldflags, opencl_lib_dir, libopencl );


%% CLODE
mex('clODEmex.cpp',[clode_path,'OpenCLResource.cpp'],[clode_path,'CLODE.cpp'],...
    debugchar,verbosechar,compflags,...
    ['-DCLODE_ROOT=\"' clode_path '\"'],...
    ['-I' spdlog_path], ['-I' clode_path], ['-I' opencl_include_dir],...
    ldflags, opencl_lib_dir, libopencl );

%% CLODEfeatures
mex('clODEfeaturesmex.cpp',[clode_path,'OpenCLResource.cpp'],[clode_path,'CLODE.cpp'],...
    [clode_path,'CLODEfeatures.cpp'],...
    debugchar,verbosechar,compflags,...
    ['-DCLODE_ROOT=\"' clode_path '\"'],...
    ['-I' spdlog_path], ['-I' clode_path], ['-I' opencl_include_dir],...
    ldflags, opencl_lib_dir, libopencl );

%% CLODEtrajectory
mex('clODEtrajectorymex.cpp',[clode_path,'OpenCLResource.cpp'],[clode_path,'CLODE.cpp'],...
    [clode_path,'CLODEtrajectory.cpp'],...
    debugchar,verbosechar,compflags,...
    ['-DCLODE_ROOT=\"' clode_path '\"'],...
    ['-I' spdlog_path], ['-I' clode_path], ['-I' opencl_include_dir],...
    ldflags, opencl_lib_dir, libopencl );