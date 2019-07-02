%test basic clODE
clear

clSinglePrecision=true;
odefile='lactotroph.ode';
[~,prob]=ode2cl(odefile,[],clSinglePrecision);

stepper='dorpri5';
% stepper='bs23'; 
% stepper='rk4'; 
vendor='nvidia';
% vendor='intel';
devicetype='default';
observer='nhood2';

sp.dt=0.1;
sp.dtmax=1.00;
sp.abstol=1e-7;
sp.reltol=1e-7;
sp.max_steps=100000;
sp.max_store=20000; %making this bigger than adaptive solver needs slows things down
sp.nout=1;

spt=sp;
% spt.dt=10;
% spt.dtmax=150.00;
% spt.abstol=1e-7;
% spt.reltol=1e-5;
% spt.max_steps=1000000;
% spt.max_store=20000; %making this bigger than adaptive solver needs slows things down
% spt.nout=1;
spt.nout=1;

%allocated number of timepoints: min( (tf-t0)/(dt*nout)+1 , max_store).
%if not enough, tf is not reached. may be too many for adaptive stepper


op.eVarIx=4;
op.fVarIx=1;
op.maxEventCount=5000;
op.minXamp=0;
op.minIMI=0;
op.nHoodRadius=0.01;
op.xUpThresh=0.3;
op.xDownThresh=0.2;
op.dxUpThresh=0;
op.dxDownThresh=0;
op.eps_dx=1e-10;

clo=clODEfeatures(prob, stepper, observer, clSinglePrecision,vendor,devicetype);

cloTraj=clODEtrajectory(prob, stepper, clSinglePrecision,vendor,devicetype);

%%
tspan=[0,10000];
nGrid=[32,32];

nPts=prod(nGrid);

mySeed=1;

p=prob.p0;

plb=[prob.par.lb];
pub=[prob.par.ub];
p1ix=find(prob.parNames=="gbk");
p2ix=find(prob.parNames=="taubk");

p1=linspace(plb(p1ix),pub(p1ix),nGrid(1));
p2=linspace(plb(p2ix),pub(p2ix),nGrid(2));
[P1,P2]=meshgrid(p1,p2);

x0=prob.x0;
X0=repmat(x0,nPts,1);

% x0lb=[prob.var.lb];
% x0ub=[prob.var.ub];
% x0lb(1)=-70;
% x0ub(1)=0;
% X0(:,1)=x0lb(1)+rand(nPts,1)*(x0ub(1)-x0lb(1));

P=repmat(p,nPts,1);
P(:,p1ix)=P1(:);
P(:,p2ix)=P2(:);

%
clo.initialize(tspan, X0(:), P(:), sp, op);
clo.seedRNG(0)

X0t=repmat(x0,1,1);
Pt=repmat(p,1,1);
cloTraj.initialize(tspan, X0t(:), Pt(:), spt);
cloTraj.tscale=1/1000;
%%
tic
clo.transient();
toc

tic
clo.transient();
toc
% 
% tic
% clo.transient();
% toc

% tspan=[0,50000];
% clo.settspan(tspan);
%%
tic
clo.features(1);
toc

tic
clo.features();
toc


%% plot

F=clo.getF();

fix=6;
f=reshape(F(:,fix),nGrid);

figure(1); clf
hi=imagesc(p1,p2,f);
hi.HitTest='off';
set(gca,'ydir','normal');
hcb=colorbar('northoutside');
xlabel(prob.parNames(p1ix));
ylabel(prob.parNames(p2ix));

title(hcb,clo.fNames{fix})

ax=gca;
axis square

nClick=3;

ax.ButtonDownFcn={@clickTrajectory,cloTraj,prob,p,x0,tspan,[p1ix,p2ix],1,2,nClick};