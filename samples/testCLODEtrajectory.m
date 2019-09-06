%test basic clODE
clear

clSinglePrecision=true;
odefile='~/Dropbox/CurrentWork/xppSrc/BJ_04a_nozzle.ode';
% odefile='C:\Users\fletcherpa\Dropbox\CurrentWork\xppSrc\BJ_04a_nozzle.ode';
% odefile='lactotroph.ode';
[~,prob]=ode2cl(odefile,[],clSinglePrecision);

stepper='dorpri5';
vendor='nvidia';
devicetype='default';

sp.dt=10;
sp.dtmax=150.00;
sp.abstol=1e-6;
sp.reltol=1e-3;
sp.max_steps=10000000;
sp.max_store=20000; %making this bigger than adaptive solver needs slows things down
sp.nout=1;

%allocated number of timepoints: min( (tf-t0)/(dt*nout)+1 , max_store).
%if not enough, tf is not reached. may be too many for adaptive stepper

tspan=[0,300000];

mySeed=1;

nPts=32;

p=prob.p0;
plb=[prob.par.lb];
pub=[prob.par.ub];

% plb=[0,5,0];
% pub=[2.5,15,2];

p(5)=3;
p(6)=.5;

% x0=[0,0,0,0];
x0=prob.x0;

p1ix=1;
p1=linspace(plb(p1ix),pub(p1ix),nPts);

X0=repmat(x0,nPts,1);
P=repmat(p,nPts,1);
P(:,p1ix)=p1(:);

clo=clODEtrajectory(prob, stepper, clSinglePrecision,vendor,devicetype);
clo.initialize(tspan, X0(:), P(:), sp);
% clot.seedRNG(0)
% tic
% clo.transient();
% toc
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
ylabel(prob.varNames(vix))
title([prob.parNames{p1ix} '=' num2str(p1(tix))])
% plot(tt,xx(:,1),'-o')
% xlim(tspan)

% figure(3)
% line(tt(2:end),diff(tt))

%%
% N=10;
% tic
% for i=1:N
% % clot.transient();
% clot.trajectory();
% end
% time=toc;
% 
% trep=time/N