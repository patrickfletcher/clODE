function plot_twopar(clo, Grid, feature, click_pars)

fig = gcf;
Grid.dim=fliplr(Grid.dim); %need to solve row-column convert
F=reshape(feature.fun(clo.getF()), Grid.dim); 

hi=imagesc(Grid.x, Grid.y, F);
hi.HitTest='off';
ax2par=gca;
ax2par.YDir='normal';
hcb=colorbar('northoutside');
xlabel(Grid.xname);
ylabel(Grid.yname);
title(hcb,feature.name)
axis square

fig.KeyPressFcn={@gridKeyPress, clo, hi, Grid, feature};

trajFig=[];

% set up a clODEtrajectory object for clickTrajectory
if exist('click_pars','var')
    nClick=click_pars.nClick;
    vars=click_pars.vars;
    
    if ~isfield(click_pars,'fig') 
        trajFig=figure('Visible',false);
    else
        trajFig=click_pars.fig;
    end
    if ~isnumeric(trajFig)
        trajFig=trajFig.Number;
    end

    trajDevice = find({clo.devices(:).type}=="CPU"); %auto selects CPU
    tspanT=clo.tspan;
    spt=clo.sp; %copy feature solver's parameters
    precision = clo.precision;
    stepper = clo.stepper;

    cloTraj=clODEtrajectory(clo.prob,precision,trajDevice,stepper);
    cloTraj.tscale=clo.tscale;
    cloTraj.tunits=clo.tunits;

    X0t=repmat(clo.prob.x0,nClick,1);
    Pt=repmat(clo.prob.p0,nClick,1);
    cloTraj.initialize(tspanT, X0t, Pt, spt);

    %attach the "clicker" to the imagesc object in Figure 1. 
    ax2par.ButtonDownFcn={@clickTrajectory,cloTraj,Grid,vars,trajFig,nClick};
end
