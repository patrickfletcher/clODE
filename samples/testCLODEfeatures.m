%test basic clODE
clear

odefile='lactotroph.ode';
precision='single';
% precision='double';
stepper='dopri5';
% stepper='bs23'; 
% stepper='rk4';

openclDevices=queryOpenCL(); %inspect this struct to see properties of OpenCL devices found

% selectedDevice=selectDevice(); %TODO: function to get a device with
% specific properties: e.g. vendor=nvidia, type=cpu, fastest clock, most
% compute units...

%Device to use for feature grid computation
%The following uses the default device: first gpu found. It also parses the
%ODEfile and writes the OpenCL code for the ODE system. 
selectedDevice=1; %autoselect: gpu>cpu
clo=clODEfeatures(odefile,precision,selectedDevice,stepper);
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

%%
tspan=[0,30000];
nGrid=[64,64];


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


clo.initialize(tspan, X0, P, sp, op);
clo.seedRNG(0)

%%

nTrans=1;
for i=1:nTrans
tic
clo.transient();
clo.shiftX0(); %sets X0 to continue from the end of the transient
toc
end

%% 
nCont=1;
for i=1:nCont
tic
clo.features();
clo.shiftX0();
toc
end

%% plot
%display list of features recorded:
% clo.fNames

%build a feature-selection function, Ffun. The following simply extracts
%feature sixth index:
fix=1; 
% fix=find(clo.fNames==fname);
fscale=1; %in case want to change feature's units
Ffun=@(F)F(:,fix)*fscale; 
ftitle=clo.fNames{fix}; %grab the feature name from the object

%another example that could be used to detect compound bursting using nhood2:
% ftitle='relative deviation of period (% max period)';
% Ffun=@(F) (F(:,2)-F(:,3))./F(:,2)*100;

F=reshape(Ffun(clo.getF()),fliplr(nGrid));

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

% initBounds=[plb(p1ix),pub(p1ix),plb(p2ix),pub(p2ix)];
% % hf.KeyPressFcn={@gridKeyPress,clo, hi, Ffun, ftitle, nGrid,[p1ix,p2ix],initBounds};
% hf.KeyPressFcn={@gridKeyPress,clo, hi, Grid, feature};

%% set up a clODEtrajectory object for clickTrajectory
% also provides keypresses to interact with the trajectories:
% g - Go: integrate over tspan from present X0
% l - Last: set X0 to Xf from previous integration, same tspan
% s - Shift: shift both X0 and t0, then integrate
% c - Continue: same as shift, but append the result to previous integration
% r - Randomize X0 within x0lb=[clo.prob.var.lb]; and x0ub=[clo.prob.var.ub];

% nClick=3; %number of mouse clicks to select parameters for trajectories
% % tspanT=tspan;
% tspanT=[0,10000];
% 
% %variables to plot (use the name)
% vars={'v'}; %one variable
% % vars={'v','c'}; %two variables
% 
% selectedDevice=1; %choose a CPU if available (device with highest clockrate)
% % clo.prob contains all the info parsed from the ODEfile above; reuse it:
% cloTraj=clODEtrajectory(clo.prob,precision,selectedDevice,stepper); 
% 
% spt=clo.sp; %copy feature solver's parameters
% % spt=clODE.solverParams(); %default params struct
% % spt.dt=0.01; %starting dt. 
% % spt.dtmax=1000.00;
% % spt.abstol=1e-7;
% % spt.reltol=1e-5;
% spt.max_steps=10000000;
% spt.max_store=20000; %number of time points allocated for trajectory result
% spt.nout=1;
% 
% %Note: large max_store will slow things down. If final time in tspan is not
% %reached, increase max_store. To query actual number of steps taken:
% % cloTraj.getNstored
% 
% X0t=repmat(x0,1,1);
% Pt=repmat(p,1,1);
% cloTraj.initialize(tspanT, X0t(:), Pt(:), spt);
% cloTraj.tscale=1/1000;
% cloTraj.tunits='s';
% 
% %attach the "clicker" to the imagesc object in Figure 1. 
% trajFig=2; %will plot trajectories in Figure 2
% ax.ButtonDownFcn={@clickTrajectory,cloTraj,p,x0,tspanT,[p1ix,p2ix],vars,trajFig,nClick};