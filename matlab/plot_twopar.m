function plot_twopar(fig, clo, Grid, feature, click_pars)
if ~isfield(Grid,'xlabel')
    Grid.xlabel=Grid.xname;
end
if ~isfield(Grid,'ylabel')
    Grid.ylabel=Grid.yname;
end

Grid.dim=fliplr(Grid.dim); %need to solve row-column convert
F=reshape(feature.fun(clo.getF()), Grid.dim); 

F(F<-1e10|F>1e10)=nan;

figure(fig)
clf reset

ax2par=axes();
hi=imagesc(ax2par, Grid.x, Grid.y, F);
hi.HitTest='off';
ax2par.YDir='normal';
hcb=colorbar('northoutside');
xlabel(Grid.xlabel);
ylabel(Grid.ylabel);
title(hcb,feature.name,'Interpreter','none')
axis square

fig.KeyPressFcn={@gridKeyPress,clo, hi, Grid, feature};

% set up a clODEtrajectory object for clickTrajectory
if exist('click_pars','var')
%     nClick=click_pars.nClick;
%     vars=click_pars.vars;
    
    X0grid=clo.getXf;
    Pgrid=clo.P(:,[Grid.xix,Grid.yix]);
    
    %set up a figure (represent by number in case we need to make it again)
    if ~isfield(click_pars,'fig') 
        click_pars.fig=figure('Visible',false);
    end
    if class(click_pars.fig)=="matlab.ui.Figure"
        click_pars.fig=click_pars.fig.Number;
    end

    trajDevice = find({clo.devices(:).type}=="CPU"); %auto selects CPU
    tspanT=clo.tspan;
    spt=clo.sp; %copy feature solver's parameters
    precision = clo.precision;
    stepper = clo.stepper;
%     stepper = 'dopri5';

    cloTraj=clODEtrajectory(clo.prob,precision,trajDevice,stepper);
    cloTraj.tscale=clo.tscale;
    cloTraj.tunits=clo.tunits;
    cloTraj.buildCL();

    X0t=repmat(clo.prob.x0,1,1);
    Pt=repmat(clo.prob.p0,1,1);
    cloTraj.initialize(tspanT, X0t(:), Pt(:), spt);

    %attach the "clicker" to the imagesc object
    ax2par.ButtonDownFcn={@clickTrajectory,cloTraj,Grid,click_pars,'k',X0grid,Pgrid};
end