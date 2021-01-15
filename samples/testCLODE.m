%test basic clODE
clear

odefile='lactotroph.ode';
precision='single';
clo=clODE(odefile,precision);
% clo.stepper='seuler'; %default='dopri5'
% clo.selectDevice(); %{'type','gpu'}, {'vendor','nvidia'}, {platID,devID}, 'maxComputeUnits','maxClock' 

%select ode problem (clo.prob), stepper, precision before building on
%current device
clo.buildCL();

%solver parameters
sp=clODE.defaultSolverParams();%create required ODE solver parameter struct
% sp.dt=0.1;
% sp.dtmax=100.00;
% sp.abstol=1e-6;
% sp.reltol=1e-3;
% sp.max_steps=10000000;

tspan=[0,1000];

%set up parameters and initial conditions
nPts=32;
p=clo.prob.p0;
x0=[0,0,0,0];

X0=repmat(x0,nPts,1);
P=repmat(p,nPts,1);


clo.initialize(tspan, X0, P, sp);
clo.seedRNG(42)
clo.transient(); %warm up the GPU for timing

%run the simulation
tic
clo.transient();
toc

%read data back from the GPU
tic
nextTspan=clo.getTspan;
xf=clo.getXf;
toc
