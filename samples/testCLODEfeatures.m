%test basic clODE
clear

openclDevices=queryOpenCL(); %inspect this struct to see properties of OpenCL devices found

odefile='lactotroph.ode';
tic
clo=clODEfeatures(odefile);
toc

clo.precision='single';
% clo.precision='double';
clo.stepper='dopri5';
% clo.stepper='bs23'; 
% clo.stepper='rk4';

%Device to use for feature grid computation
%The following uses the default device: first gpu found. It also parses the
%ODEfile and writes the OpenCL code for the ODE system. 
% selectedDevice=[]; %autoselect: gpu>cpu
% tic
% clo=clODEfeatures(odefile,precision,selectedDevice,stepper);
% toc
clo.getProgramString

tic
clo.buildCL();
toc

clo.printStatus

%set properties 
% clo.stepper='rk4'; %default='dorpri5'

%solver parameters`
sp=clODE.defaultSolverParams();%create required ODE solver parameter struct
sp.dt=0.5;
sp.dtmax=100;
sp.abstol=1e-6;
sp.reltol=1e-4; %nhood2 may require fairly strict reltol
sp.max_steps=1e7;

op=clODEfeatures.defaultObserverParams(); %create required observer parameter struct
op.maxEventCount=10000; %stops if this many events found {localmax, nhood2}
op.eps_dx=0e-7; %for checking for min/max
op.minXamp=0.00; %don't count event if global (max x - min x) is too small
op.minIMI=0.00; %don't count event if global (max x - min x) is too small

% clo.observer='basic'; %records the extent (max/min) of a variable and its slope
% clo.observer='basicall'; %same as above but for all variables

% clo.observer='localmax'; %features derived from local maxima and minima only
% fname="max IMI";

% clo.observer='nhood1'; %start point is first local min of variable eVarIx
% op.eVarIx=4; %variable used for deciding centerpoint of neighborhood
% op.fVarIx=1; %feature detection variable
% op.nHoodRadius=.2; %size of neighborhood

% clo.observer='nhood2'; %"Poincare ball": period detection by trajectory returning to within a neighborhood of a specific point in state space
% op.eVarIx=4; %nhood2: variable used for deciding centerpoint of neighborhood
% op.fVarIx=1; %feature detection variable
% op.nHoodRadius=.2; %size of neighborhood {nhood2} 
% op.xDownThresh=0.05; %selecting neighborhood centerpoint: first time eVarIx drops below this fraction of its amplitude 

% clo.observer='thresh2'; %event detection and features both measured in variable fVarIx
% op.fVarIx=1;
% %for constructing up/down thresholds:
% op.xUpThresh=0.3; %must provide xUpThresh at least
% op.xDownThresh=0.1; %xDownThresh=0 => use same as xUpThresh
% op.dxUpThresh=0.; %dxUpThresh=0 => don't use
% op.dxDownThresh=0.; %dxUpThresh=0 => use same as dxUpThresh
% % % fname = "mean peaks";

%display list of features that will be computed:
% clo.fNames

%%
tspan=[0,30000];
nGrid=[64,64];
nPts=prod(nGrid);

p1ix=find(clo.prob.parNames=="gbk");
p2ix=find(clo.prob.parNames=="taubk");

p=clo.prob.p0; plb=[clo.prob.par.lb]; pub=[clo.prob.par.ub];
p1=linspace(plb(p1ix),pub(p1ix),nGrid(1));
p2=linspace(plb(p2ix),pub(p2ix),nGrid(2));
[P1,P2]=meshgrid(p1,p2);

P=repmat(p,nPts,1);
P(:,p1ix)=P1(:);
P(:,p2ix)=P2(:);

x0lb=[clo.prob.var.lb];
x0ub=[clo.prob.var.ub];
X0=generate_pointset('random',x0lb,x0ub,nPts);

clo.initialize(tspan, X0, P, sp, op);

%% first run
%run a transient
tic
clo.transient();
clo.shiftX0(); %sets X0 to continue from the end of the transient
toc

%compute features
tic
clo.features(1);
toc

% get data from the device
F=clo.getF();

% plot
fix=1; %feature to plot
ff=reshape(F(:,fix),fliplr(nGrid));

hf=figure(1); clf
hi=imagesc(p1,p2,ff);
hi.HitTest='off';
ax=gca;
ax.YDir='normal';
hcb=colorbar('northoutside');
xlabel(clo.prob.parNames(p1ix));
ylabel(clo.prob.parNames(p2ix));
title(hcb,clo.fNames{fix})
axis square
