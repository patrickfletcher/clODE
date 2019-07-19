function clickTrajectory(clickAxis,~,clo,ivp,pin,y0in,tspan,parIx,varIx,tracjectoryFigID,nClick,markerColor)
%A mouse click inside the axes will trigger a parameter
%selection crosshair. A second click selects a parameter
%combination to simulate, and resulting simulation is displayed
%in a new figure window.

%could find closest solved instance, get y0 for beginning of
%feature detection interval

if ~exist('varIx','var'), varIx=1; end
if ~exist('nClick','var'), nClick=1; end
if nClick>6, error('I refuse to do more than 6 clicks'); end

if ~exist('markerColor','var'), markerColor='k';end

% y0=repmat([ivp.var(:).value],nClick,1);
% p=repmat([ivp.par(:).value],nClick,1);
y0=repmat(y0in(:)',nClick,1);
p=repmat(pin(:)',nClick,1);
vNames={ivp.var(:).name};
pNames={ivp.par(:).name};

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

clo.setProblemData(y0(:),p(:));
clo.settspan(tspan);
% clo.transient();

if ~exist('tracjectoryFigID','var'), tracjectoryFigID=figure();end
hf=figure(tracjectoryFigID);

hf.KeyPressFcn=@keypress;

clf
[X,T]=integrate();
plotTrajectories(X,T);

    function [X,T,nStored]=integrate()
        clo.trajectory();
        xx=clo.getX();
        %         AUX=clo.getAux();
        tt=clo.getT()*clo.tscale;
        nStored=clo.getNstored();
        for i=1:nClick
            T{i}=tt(1:nStored(i),i);
            X{i}=xx(1:nStored(i),:,i);
        end
    end

    function keypress(src,evt)
        switch(evt.Key)
            case 'g'
                [X,T]=integrate();
                plotTrajectories(X,T);
                
            case 'c'
                [xx,tt]=integrate();
                for i=1:nClick
                    T{i}=[T{i};tt{i}+T{i}(end)];
                    X{i}=[X{i};xx{i}];
                end
                plotTrajectories(X,T);
        end
    end

    function plotTrajectories(X,T)
        
        for i=1:nClick
            
            t=T{i};
            x=X{i};
            
            %Time plot
            ax=subplot(nClick,1,i);
            if numel(varIx)==1
                plot(t,x(:,varIx),'k');
                xlim([t(1),t(end)]);
                ylabel(vNames{varIx});
                YLIM=ylim();
                xmax=max(x(:,varIx(1)));
                xmin=min(x(:,varIx(1)));
                ylim([min(YLIM(1),xmin-0.01*abs(xmin)), max(YLIM(2),xmax+0.01*abs(xmax))]);
                
            elseif numel(varIx)==2
                %Time plot
                yyaxis left
                plot(t,x(:,varIx(1)));
                xlim([t(1),t(end)]);
                YLIM=ylim();
                xmax=max(x(:,varIx(1)));
                xmin=min(x(:,varIx(1)));
                ylim([min(YLIM(1),xmin-0.01*abs(xmin)), max(YLIM(2),xmax+0.01*abs(xmax))]);
                
                ylabel(vNames{varIx(1)});
                
                yyaxis right
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
            title([pNames{parIx(1)} '=' num2str(p(i,parIx(1)),2) ', ' pNames{parIx(2)} '=' num2str(p(i,parIx(2)),2)]);
        end
    end
end