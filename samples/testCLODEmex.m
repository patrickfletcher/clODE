%test basic clODEmex

clSinglePrecision=false;
[rhsfile,prob]=ode2cl('lactotroph.ode',[],clSinglePrecision);

vendor=0; %any
devicetype=0; %default
stepper=2; %2-rk4

sp.dt=0.5;
sp.dtmax=1.00;
sp.abstol=0.01;
sp.reltol=0.01;
sp.max_steps=10000000;
sp.max_store=10000;
sp.nout=50;

tspan=[0,500];

mySeed=1;

nPts=4096;

p=[1,5,1];
x0=[0,0,0,0];

X0=repmat(x0,nPts,1);
P=repmat(p,nPts,1);

clo=clODEmex('new',prob, stepper, clSinglePrecision,vendor, devicetype);
clODEmex('initialize',clo, tspan, X0(:), P(:), sp);
clODEmex('seedrng',clo,0);

tic
clODEmex('transient',clo);
toc


clODEmex('settspan',clo, tspan);
clODEmex('setproblemdata',clo, X0(:), P(:));
tic
clODEmex('transient',clo);
toc


tspan=clODEmex('getTspan',clo);
xf=clODEmex('getX0',clo);

XF=reshape(xf,nPts,prob.nVar);
clODEmex('delete',clo);