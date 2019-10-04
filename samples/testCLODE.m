%test basic clODE
clear

odefile='lactotroph.ode';

stepper='dorpri5';
vendor='nvidia';
devicetype='default';
clSinglePrecision=true;

[rhsfile,prob]=ode2cl(odefile,[],clSinglePrecision);

sp.dt=0.1;
sp.dtmax=100.00;
sp.abstol=1e-6;
sp.reltol=1e-3;
sp.max_steps=10000000;
sp.max_store=10000;
sp.nout=50;

tspan=[0,1000];

nPts=32;
p=prob.p0;
x0=[0,0,0,0];

X0=repmat(x0,nPts,1);
P=repmat(p,nPts,1);

clo=clODE(prob, stepper, clSinglePrecision,vendor,devicetype);
clo.initialize(tspan, X0(:), P(:), sp);
clo.seedRNG(42)
% clo.transient();

tic
clo.transient();
toc

tic
nextTspan=clo.getTspan;
xf=clo.getXf;
auxf=clo.getAuxf;
toc
