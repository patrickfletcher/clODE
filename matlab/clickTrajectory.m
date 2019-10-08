function clickTrajectory(clickAxis,~,clo,pin,x0in,tspan,parIx,vars,tracjectoryFigID,nClick,markerColor)
%A mouse click inside the axes will trigger a parameter
%selection crosshair. A second click selects a parameter
%combination to simulate, and resulting simulation is displayed
%in a new figure window.

%could find closest solved instance, get x0 for beginning of
%feature detection interval

if ~exist('vars','var'), vars=clo.prob.varNames(1); end
if ischar(vars), vars={vars}; end
if ~exist('nClick','var'), nClick=1; end
if nClick>6, error('I refuse to do more than 6 clicks'); end

if ~exist('markerColor','var'), markerColor='k';end

x0=repmat(x0in(:)',nClick,1);
p=repmat(pin(:)',nClick,1);

varIsAux=false(size(vars));
for v=1:length(vars) 
    tmpIx=find(clo.prob.varNames==string(vars{v})); 
    if isempty(tmpIx)
        tmpIx=find(clo.prob.auxNames==string(vars{v}));
        if isempty(tmpIx)
            error(['unknown variable: ' vars(v)]);
        else
            varIsAux(v)=true;
            varIx(v)=tmpIx;
        end 
    else
        varIx(v)=tmpIx;
    end
end

markers='osdp^v';

oldMarks=findobj('tag','marks');
if ~isempty(oldMarks)
    delete(oldMarks);
end

for n=1:nClick
    [px,py]=ginput(1);
    p(n,parIx(1))=px;
    p(n,parIx(2))=py;
    h(n)=line(px,py,'marker',markers(n),'linestyle','none','color',markerColor,'linewidth',1,'tag','marks');
end

clo.setProblemData(x0(:),p(:));
clo.settspan(tspan);
% clo.transient();

if ~exist('tracjectoryFigID','var'), tracjectoryFigID=figure();end
hf=figure(tracjectoryFigID);

hf.KeyPressFcn=@keypress;

clf
[X,T,AUX]=integrate();
plotTrajectories(X,T,AUX);

    function [X,T,AUX]=integrate()
        clo.trajectory();
        xx=clo.getX();
        aux=clo.getAux();
        tt=clo.getT()*clo.tscale;
        nStored=clo.getNstored();
        for i=1:nClick
            T{i}=tt(1:nStored(i),i);
            X{i}=xx(1:nStored(i),:,i);
            AUX{i}=aux(1:nStored(i),:,i);
        end
    end

    function keypress(src,evt)
        switch(evt.Key)
            case 'g'
                clo.settspan(tspan);
                [X,T,AUX]=integrate();
                plotTrajectories(X,T,AUX);
                
            case 'c'
                clo.shiftX0();
                [xx,tt,aux]=integrate();
                for i=1:nClick
                    T{i}=[T{i};tt{i}+T{i}(end)];
                    X{i}=[X{i};xx{i}];
                    AUX{i}=[AUX{i};aux{i}];
                end
                plotTrajectories(X,T,AUX);
                
            case 'r' %randomize ICs
                clo.settspan(tspan);
                x0lb=[clo.prob.var.lb];
                x0ub=[clo.prob.var.ub];
                x0=x0lb+rand(nClick,length(x0lb)).*(x0ub-x0lb);
                clo.setX0(x0(:));
                [X,T,AUX]=integrate();
                plotTrajectories(X,T,AUX);
        end
    end

    function plotTrajectories(X,T,AUX)
        
        
        for i=1:nClick
            
            t=T{i};
            
            %Time plot
            ax=subplot(nClick,1,i);
            if numel(varIx)==1
                
                if varIsAux
                    x=AUX{i};
                    vNames=clo.prob.auxNames;
                else
                    x=X{i};
                    vNames=clo.prob.varNames;
                end
                
                plot(t,x(:,varIx),'k');
                xlim([t(1),t(end)]);
                ylabel(vNames{varIx});
                YLIM=ylim();
                xmax=max(x(:,varIx));
                xmin=min(x(:,varIx));
                ylim([min(YLIM(1),xmin-0.01*abs(xmin)), max(YLIM(2),xmax+0.01*abs(xmax))]);
                
            elseif numel(varIx)==2
                %Time plot
                yyaxis left
                if varIsAux(1)
                    x=AUX{i};
                    vNames=clo.prob.auxNames;
                else
                    x=X{i};
                    vNames=clo.prob.varNames;
                end
                plot(t,x(:,varIx(1)));
                xlim([t(1),t(end)]);
                YLIM=ylim();
                xmax=max(x(:,varIx(1)));
                xmin=min(x(:,varIx(1)));
                ylim([min(YLIM(1),xmin-0.01*abs(xmin)), max(YLIM(2),xmax+0.01*abs(xmax))]);
                
                ylabel(vNames{varIx(1)});
                
                yyaxis right
                if varIsAux(2)
                    x=AUX{i};
                    vNames=clo.prob.auxNames;
                else
                    x=X{i};
                    vNames=clo.prob.varNames;
                end
                plot(t,x(:,varIx(2)));
                xlim([t(1),t(end)]);
                YLIM=ylim();
                xmax=max(x(:,varIx(2)));
                xmin=min(x(:,varIx(2)));
                ylim([min(YLIM(1),xmin-0.01*abs(xmin)), max(YLIM(2),xmax+0.01*abs(xmax))]);
                ylabel(vNames{varIx(2)});
            end
            
            box off
            
            if i<nClick
                ax.XTickLabel=[];
            else
                xlabel(['t' ' (' clo.tunits ')']);
            end
            
            XLIM=xlim();
            
            line((XLIM(1)+XLIM(end))*0.02,YLIM(2)*0.9,'marker',markers(i),'linestyle','none','color',markerColor,'linewidth',1,'tag','marks');
            title([clo.prob.parNames{parIx(1)} '=' num2str(p(i,parIx(1)),2) ', ' clo.prob.parNames{parIx(2)} '=' num2str(p(i,parIx(2)),2)]);
        end
    end
end