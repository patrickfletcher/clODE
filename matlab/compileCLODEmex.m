%test for new mex version of clODE:
clear

% %%%%%%%%%%%%%%%% CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ismac
    % Code to run on Mac plaform
elseif isunix
    % Code to run on Linux plaform
    opencl_include_dir = '/usr/local/cuda-7.5/targets/x86_64-linux/include';
    opencl_lib_dir = '/usr/lib64/nvidia';
elseif ispc
    % Code to run on Windows platform 
    opencl_include_dir = 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.1/include';
    opencl_lib_dir = 'C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.1/lib/x64';
else
    disp('Cannot recognize platform')
end

clode_path=pwd; 
clode_path=clode_path(1:find(clode_path==filesep,1,'last'));
clode_path=[clode_path,'src/'];
clode_path=strrep(clode_path,'\','/'); %stupid windows backslash filesep

debugchar='';
% debugchar='-g';

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mex('queryOpenCL.cpp',[clode_path,'OpenCLResource.cpp'],...
    debugchar,...
    ['-I' clode_path], ['-I' opencl_include_dir], ['-L' opencl_lib_dir], '-lOpenCL' );


%% CLODE
mex('clODEmex.cpp',[clode_path,'OpenCLResource.cpp'],[clode_path,'CLODE.cpp'],...
    debugchar,...
    ['-DCLODE_ROOT=\"' clode_path '\"'],...
    ['-I' clode_path], ['-I' opencl_include_dir], ['-L' opencl_lib_dir], '-lOpenCL' );

%%

mex('clODEtrajectorymex.cpp',[clode_path,'OpenCLResource.cpp'],[clode_path,'CLODE.cpp'],...
    [clode_path,'CLODEtrajectory.cpp'],...
    debugchar,...
    ['-DCLODE_ROOT=\"' clode_path '\"'],...
    ['-I' clode_path], ['-I' opencl_include_dir], ['-L' opencl_lib_dir], '-lOpenCL' );

%%

mex('clODEfeaturesmex.cpp',[clode_path,'OpenCLResource.cpp'],[clode_path,'CLODE.cpp'],...
    [clode_path,'CLODEfeatures.cpp'],...
    debugchar,...
    ['-DCLODE_ROOT=\"' clode_path '\"'],...
    ['-I' clode_path], ['-I' opencl_include_dir], ['-L' opencl_lib_dir], '-lOpenCL' );