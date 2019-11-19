%test basic clODE
clear

odefile='lactotroph.ode';
precision='single';
clo=clODEtrajectory(odefile,precision);
% clo.stepper='rk4'; %default='dopri5'

%solver parameters
sp=clODE.solverParams();%create required ODE solver parameter struct
% sp.dt=0.1;
% sp.dtmax=100.00;
% sp.abstol=1e-6;
% sp.reltol=1e-3;
% sp.max_steps=10000000;

tspan=[0,3000];

%set up parameters and initial conditions
nPts=32;

X0=repmat(clo.prob.x0,nPts,1);
P=repmat(clo.prob.p0,nPts,1);


clo.initialize(tspan, X0, P, sp);
clo.seedRNG(42)
clo.transient(); %warm up the GPU for timing
%%
tic
clo.trajectory();
toc

tic
t=clo.getT();
x=clo.getX();
dx=clo.getDx();
aux=clo.getAux();
nSteps=clo.getNsteps();
nStored=clo.getNstored();
toc


%%
figure(1)
tix=32;
vix=1;
thisNstore=nStored(tix)
tt=t(1:thisNstore,tix);
xx=x(1:thisNstore,vix,tix);
auxx=aux(1:thisNstore,tix);
plot(tt,xx(:,1))
xlabel('t')
ylabel(clo.prob.varNames(vix))

%% show dt(t)
figure(2)
tix=32;
vix=1;
tt=t(1:thisNstore,tix);
plot(tt(1:end-1),diff(tt))
xlabel('t')
ylabel('dt')