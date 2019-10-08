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
clo=clODEfeatures(odefile,precision); 

%set properties 
% clo.stepper='rk4'; %default='dorpri5'
clo.observer='localmax';
% clo.observer='nhood2';

%device to use for trajectory computation
% clo.prob contains all the info parsed from the ODEfile above. We can
% reuse that info by passing it in as first argument.
selectedDevice=2; %select a specific device
cloTraj=clODEtrajectory(clo.prob,precision,selectedDevice); 



sp.dt=0.01;
sp.dtmax=100.00;
sp.abstol=1e-7;
sp.reltol=1e-5;
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

op.eVarIx=1;
op.fVarIx=1;
op.maxEventCount=5000;
op.minXamp=1;
op.minIMI=0; %not implemented: 10*sp.dtmax
op.nHoodRadius=0.1;
op.xUpThresh=0.3;
op.xDownThresh=0.05;
op.dxUpThresh=0;
op.dxDownThresh=0;
op.eps_dx=1e-7;

%%
tspan=[0,1000];
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

X0t=repmat(x0,1,1);
Pt=repmat(p,1,1);
cloTraj.initialize(tspan, X0t(:), Pt(:), spt);
cloTraj.tscale=1/1000;
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

F=clo.getF();

fix=1;
f=reshape(F(:,fix),nGrid);

figure(1); clf
hi=imagesc(p1,p2,f);
hi.HitTest='off';
ax=gca;
ax.YDir='normal';
hcb=colorbar('northoutside');
xlabel(clo.prob.parNames(p1ix));
ylabel(clo.prob.parNames(p2ix));
title(hcb,clo.fNames{fix})
axis square

nClick=3;

%variables to plot (use the name)
vars={'v'}; %one variable
% vars={'v','c'}; %two variables

ax.ButtonDownFcn={@clickTrajectory,cloTraj,p,x0,tspan,[p1ix,p2ix],vars,2,nClick};