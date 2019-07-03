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
[X,T,nStored]=integrate();
plotTrajectories(X,T,nStored);

    function [X,T,nStored]=integrate()
        clo.trajectory();
        X=clo.getX();
        AUX=clo.getAux();
        T=clo.getT()*clo.tscale;
        nStored=clo.getNstored();
    end

    function keypress(src,evt)
        switch(evt.Key)
            case 'c'
                [X,T,nStored]=integrate();
                plotTrajectories(X,T,nStored);
        end 
    end

    function plotTrajectories(X,T,nStored)
        
        for i=1:nClick
            
            t=T(:,i);
            x=X(:,:,i);
            
            t=t(1:nStored(i));
            
            x=x(1:nStored(i),:);
            
            
            %Time plot
            ax=subplot(nClick,1,i);
            plot(t,x(:,varIx),'k')
            xlim([t(1),t(end)]);
            
            % To do features, need to have mRHSfun...
            %     if doFeatures
            %
            %         dy=zeros(length(t),nVar);
            %         aux=zeros(length(t),nAux);
            %         for i=1:length(t)
            %             [tmpdy, tmpaux]=mRHSfun(t(i),thisY(i,:));
            %             dy(i,:)=tmpdy(:)';
            %         end
            %
            %         thisFmatlab=observer_plateau(t,y,dy,aux,Times,obspars);
            %         Feats(:,i)=thisFmatlab(:);
            %     end
            
            if i<nClick
                ax.XTickLabel=[];
            end
            
            title([pNames{parIx(1)} '=' num2str(p(i,parIx(1))) ', ' pNames{parIx(2)} '=' num2str(p(i,parIx(2)))]);
            xlabel('t'); ylabel(vNames{varIx});
            
        end
    end

end