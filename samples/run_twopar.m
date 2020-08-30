%first script for testing input encoding to ca signals.
clear

gridDevice=1; %select a specific device: grid should use GPU (device with most "cores")
trajDevice=2; %for trajectories, choose a CPU if available (device with highest clockrate)
precision='single';
% stepper='dopri5'; %default
stepper='rk4';
Grid.dim=[32,64];

odefile='lactotroph.ode';
clo=clODEfeatures(odefile,precision,gridDevice,stepper);

% clo.getAvailableSteppers %lists available time-stepping methods

eVar="c"; %event detection variable
fVar="v"; %feature variable
traj_vars={'v','c'}; %for clicktraj
Grid.xname="gbk";
Grid.yname="taubk";
clo.tscale=1; %ms to s
clo.tunits='ms';
tspan=[0,30000];

%solver parameters`
sp=clODE.defaultSolverParams();%create required ODE solver parameter struct
sp.dt=clo.prob.opt.dt;
sp.max_store=1e6;
% sp.dtmax=100;
sp.abstol=1e-6;
sp.reltol=1e-4; %nhood2 may require fairly strict reltol
% sp.max_steps=2000;

eVarIx=find(clo.prob.varNames==eVar);
fVarIx=find(clo.prob.varNames==fVar);

%observer paramters
op=clODEfeatures.defaultObserverParams(); %create required observer parameter struct
op.maxEventCount=10000; %stops if this many events found {localmax, nhood2}
op.eps_dx=1e-7; %for checking for min/max
op.minXamp=0.01; %don't count event if global (max x - min x) is too small

% op.fVarIx=fVarIx;
% clo.observer='basic'; %records the extent (max/min) of a variable and its slope
clo.observer='basicall'; %same as above but for all variables

% clo.observer='localmax'; %features derived from local maxima and minima only
% op.fVarIx=fVarIx;

% clo.observer='nhood1'; %start point is first local min of variable eVarIx
% op.eVarIx=eVarIx; %variable used for deciding centerpoint of neighborhood
% op.fVarIx=fVarIx; %feature detection variable
% op.nHoodRadius=.2; %size of neighborhood {nhood2} 

% clo.observer='nhood2'; %"Poincare ball": period detection by trajectory returning to within a neighborhood of a specific point in state space
% op.eVarIx=eVarIx; %nhood2: variable used for deciding centerpoint of neighborhood
% op.fVarIx=fVarIx; %feature detection variable
% op.nHoodRadius=.2; %size of neighborhood {nhood2} 
% op.xDownThresh=0.5; %selecting neighborhood centerpoint: first time eVarIx drops below this fraction of its amplitude 

% clo.observer='thresh2'; %event detection and features both measured in variable fVarIx
% op.fVarIx=fVarIx;
% %for constructing up/down thresholds:
% op.xUpThresh=0.5; %must provide xUpThresh at least
% op.dxUpThresh=0.; %dxUpThresh=0 => don't use
% op.xDownThresh=0.5; %xDownThresh=0 => use same as xUpThresh
% op.dxDownThresh=0.; %dxUpThresh=0 => use same as dxUpThresh




%% 2-par grid

nPts=prod(Grid.dim);

p=clo.prob.p0;
plb=[clo.prob.par.lb];
pub=[clo.prob.par.ub];
Grid.xix = find(clo.prob.parNames==Grid.xname);
Grid.yix = find(clo.prob.parNames==Grid.yname);
Grid.x=linspace(plb(Grid.xix),pub(Grid.xix),Grid.dim(1));
Grid.y=linspace(plb(Grid.yix),pub(Grid.yix),Grid.dim(2));
[X,Y]=meshgrid(Grid.x,Grid.y);

Grid.xtype = 'par';
Grid.ytype = 'par';
P=repmat(p,nPts,1);
P(:,Grid.xix)=X(:);
P(:,Grid.yix)=Y(:);

x0=clo.prob.x0;
% X0=repmat(x0,nPts,1);
rng(1);
x0lb=[clo.prob.var.lb];
x0ub=[clo.prob.var.ub];
X0=x0lb+rand(nPts,length(x0lb)).*(x0ub-x0lb);

clo.initialize(tspan, X0, P, sp, op);


%%
nTrans=2;
for i=1:nTrans
tic
clo.transient();
clo.shiftX0(); %sets X0 to continue from the end of the transient
toc
end

tic
clo.features(1);
clo.shiftX0();
toc

nCont=1;
for i=1:nCont
tic
clo.features(0);
clo.shiftX0();
toc
end



%% plot
hf1=figure(1); clf
hf1.Units='inches';
hf1.Position(3:4)=[4,4];
hf1.PaperPosition=hf1.Position;
hf1.PaperPosition(1:2)=0;
hf1.PaperSize=hf1.Position(3:4);

fname="max v";  %clo.fNames contains list of available features
fix=find(clo.fNames == fname);
fscale=clo.tscale;
% fscale=1/1000;
% fscale=1/1000/60; %in case want to change feature's units
feature.fun=@(F)F(:,fix)*fscale; 
feature.name=clo.fNames{fix}; %grab the feature name from clo object

% feature.name='c amplitude';
% maxcix = find(clo.fNames == "max c");
% mincix = find(clo.fNames == "min c");
% feature.fun=@(F) F(:,maxcix)-F(:,mincix);

% feature.name='relative deviation of period (% max period)';
% feature.fun=@(F) (F(:,2)-F(:,3))./F(:,2)*100;
% 
% feature.name='mean upDuration / mean period';
% feature.fun=@(F) (F(:,10)./F(:,4));

click_pars.nClick=3;
click_pars.vars=traj_vars;
click_pars.fig=2;

plot_twopar(clo, Grid, feature, click_pars)

