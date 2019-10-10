%test basic clODE
clear

odefile='lactotroph.ode';
precision='single';
stepper='dorpri5';
% stepper='bs23'; 
% stepper='rk4';

openclDevices=queryOpenCL(); %inspect this struct to see properties of OpenCL devices found

% selectedDevice=selectDevice(); %TODO: function to get a device with
% specific properties: e.g. vendor=nvidia, type=cpu, fastest clock, most
% compute units...

%Device to use for feature grid computation
%The following uses the default device: first one found. It also parses the
%ODEfile and writes the OpenCL code for the ODE system. 
selectedDevice=1; %select a specific device
clo=clODEfeatures(odefile,precision,selectedDevice,stepper); 

%set properties 
% clo.stepper='rk4'; %default='dorpri5'

% clo.observer='basic'; %records the extent (max/min) of a variable and its slope
% clo.observer='basicall'; %same as above but for all variables
% clo.observer='localmax'; %features derived from local maxima and minima only
clo.observer='nhood2'; %"Poincare ball": period detection by trajectory returning to within a neighborhood of a specific point in state space


%solver parameters
sp=clODE.solverParams();%create required ODE solver parameter struct
% sp.dt=1;
% sp.dtmax=150.00;
sp.abstol=1e-9;
sp.reltol=1e-7; %nhood2 may require fairly strict reltol
% sp.max_steps=100000;

op=clODEfeatures.observerParams(); %create required observer parameter struct
op.fVarIx=1; %feature detection variable. used in: {basic, localmax, nhood2}
op.maxEventCount=1000; %stops if this many events found {localmax, nhood2}
op.minXamp=1; %don't record oscillation features if units of variable fVarIx {nhood2}
op.nHoodRadius=0.1; %size of neighborhood {nhood2} 
op.xDownThresh=0.2; %selecting neighborhood centerpoint: first time fVarIx drops below this fraction of its amplitude {nhood2}
op.eps_dx=1e-6; %for checking for min/max

%%
tspan=[0,10000];
nGrid=[32,32];

nPts=prod(nGrid);

mySeed=1;

p=clo.prob.p0;

plb=[clo.prob.par.lb];
pub=[clo.prob.par.ub];
p1ix=find(clo.prob.parNames=="gbk");
p2ix=find(clo.prob.parNames=="taubk");

p1=linspace(plb(p1ix),pub(p1ix),nGrid(1));
p2=linspace(plb(p2ix),pub(p2ix),nGrid(2));
[P1,P2]=meshgrid(p1,p2);

x0=clo.prob.x0;
X0=repmat(x0,nPts,1);

% x0lb=[prob.var.lb];
% x0ub=[prob.var.ub];
% x0lb(1)=-70;
% x0ub(1)=0;
% X0(:,1)=x0lb(1)+rand(nPts,1)*(x0ub(1)-x0lb(1));

P=repmat(p,nPts,1);
P(:,p1ix)=P1(:);
P(:,p2ix)=P2(:);


clo.initialize(tspan, X0(:), P(:), sp, op);
clo.seedRNG(0)

%%
tic
clo.transient();
toc

%%
clo.shiftX0(); %sets X0 to continue from the end of the transient

tic
clo.features();
toc

%% plot
%display list of features recorded:
% clo.fNames

%build a feature-selection function, Ffun. The following simply extracts
%feature sith index fix:
fix=7; 
fscale=1; %in case want to change feature's units
Ffun=@(F)F(:,fix)*fscale; 
ftitle=clo.fNames{fix}; %grab the feature name from the object

%another example that could be used to detect compound bursting using nhood2:
% ftitle='relative deviation of period (% max period)';
% Ffun=@(F) (F(:,2)-F(:,3))./F(:,2)*100;

F=reshape(Ffun(clo.getF()),nGrid);

hf=figure(1); clf
hi=imagesc(p1,p2,F);
hi.HitTest='off';
ax=gca;
ax.YDir='normal';
hcb=colorbar('northoutside');
xlabel(clo.prob.parNames(p1ix));
ylabel(clo.prob.parNames(p2ix));
title(hcb,ftitle)
axis square

%attach a keypress function to interact with the grid of simulations
% t - Transient: integrate forward without feature detection (will display state variable at time=tf)
% g - Go: start the feature detector and integrate over tspan
% c - Continue: extend by another tspan and continue feature detection
% r - Randomize X0 within x0lb=[clo.prob.var.lb]; and x0ub=[clo.prob.var.ub];
%
% for example, do 'r', then 't' a few times, then 'g', then 'c' until
% the heatmap stabilizes. The idea is to get to a steady-state trajectory

hf.KeyPressFcn={@gridKeyPress,clo, hi, Ffun, ftitle, nGrid};

%% set up a clODEtrajectory object for clickTrajectory
% also provides keypresses to interact with the trajectories:
% g - Go: integrate over tspan from present X0
% l - Last: set X0 to Xf from previous integration, same tspan
% s - Shift: shift both X0 and t0, then integrate
% c - Continue: same as shift, but append the result to previous integration
% r - Randomize X0 within x0lb=[clo.prob.var.lb]; and x0ub=[clo.prob.var.ub];

nClick=3; %number of mouse clicks to select parameters for trajectories

%variables to plot (use the name)
vars={'v'}; %one variable
% vars={'v','c'}; %two variables

selectedDevice=2; %choose a CPU if available (device with highest clockrate)
% clo.prob contains all the info parsed from the ODEfile above; reuse it:
cloTraj=clODEtrajectory(clo.prob,precision,selectedDevice,stepper); 

spt=clODE.solverParams(); %default params struct
% spt.dt=10;
% spt.dtmax=150.00;
% spt.abstol=1e-7;
% spt.reltol=1e-5;
% spt.max_steps=1000000;
spt.max_store=20000; %allocated number of timepoints: min( (tf-t0)/(dt*nout)+1 , sp.max_store)
spt.nout=1;


X0t=repmat(x0,1,1);
Pt=repmat(p,1,1);
cloTraj.initialize(tspan, X0t(:), Pt(:), spt);
cloTraj.tscale=1/1000;

%attach the "clicker" to the imagesc object in Figure 1. 
trajFig=2; %will plot trajectories in Figure 2
ax.ButtonDownFcn={@clickTrajectory,cloTraj,p,x0,tspan,[p1ix,p2ix],vars,trajFig,nClick};