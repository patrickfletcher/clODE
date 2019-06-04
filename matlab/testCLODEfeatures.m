%test basic clODE
clear

clSinglePrecision=true;
odefile='lactotroph.ode';
[~,prob]=ode2cl(odefile,[],clSinglePrecision);

stepper='dorpri5';
% stepper='rk4'; 
vendor='nvidia';
% vendor='intel';
devicetype='default';
observer='localmax';

sp.dt=0.1;
sp.dtmax=10.00;
sp.abstol=1e-7;
sp.reltol=1e-5;
sp.max_steps=100000;
sp.max_store=20000; %making this bigger than adaptive solver needs slows things down
sp.nout=1;

%allocated number of timepoints: min( (tf-t0)/(dt*nout)+1 , max_store).
%if not enough, tf is not reached. may be too many for adaptive stepper


op.eVarIx=1;
op.fVarIx=1;
op.maxEventCount=5000;
op.minXamp=0;
op.minIMI=0;
op.nHoodRadius=0.01;
op.xUpThresh=0.3;
op.xDownThresh=0.2;
op.dxUpThresh=0;
op.dxDownThresh=0;
op.eps_dx=1e-6;

clo=clODEfeatures(prob, stepper, observer, clSinglePrecision,vendor,devicetype);

%%
tspan=[0,100000];
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

%%
tic
clo.transient();
toc
%

% tspan=[0,50000];
% clo.settspan(tspan);

tic
clo.features();
toc

F=clo.getF();

%% plot

fix=2;
f=reshape(F(:,fix),nGrid);

figure(1); clf
imagesc(p1,p2,f)
set(gca,'ydir','normal');
hcb=colorbar('northoutside');
xlabel(prob.parNames(p1ix));
ylabel(prob.parNames(p2ix));

title(hcb,clo.fNames{fix})

axis square

drawnow