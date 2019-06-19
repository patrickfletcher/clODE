function [t,x,dY]=clickTrajectory(clickAxis,~,clo,ivp,pin,y0in,tspan,parIx,varIx,tracjectoryFigID,nClick,markerColor)
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

for i=1:nClick
    [px,py]=ginput(1);
    p(i,parIx(1))=px;
    p(i,parIx(2))=py;
    h(i)=line(px,py,'marker',markers(i),'linestyle','none','color',markerColor,'linewidth',1,'tag','marks');
end

clo.setProblemData(y0(:),p(:));
clo.settspan(tspan);
clo.transient();
clo.trajectory();
X=clo.getX;
t=linspace(tspan(1),tspan(2),length(X(:,1)));

if ~exist('tracjectoryFigID','var'), tracjectoryFigID=figure();end
figure(tracjectoryFigID)
clf
for i=1:nClick
    
    x=X(:,i:nClick:end);
    
    tix=find(t>=tspan(1),1,'first'):length(t);

    %Time plot

    subplot(nClick,1,i)
    
    plot(t(tix),x(tix,varIx),'k')
%     plot(t(tix)-t(tix(1)),x(tix,varIx),'k')  %start at t=0
    
%     ylim([sys.opt.ylo,sys.opt.yhi])
    
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

    title([pNames{parIx(1)} '=' num2str(p(i,parIx(1))) ', ' pNames{parIx(2)} '=' num2str(p(i,parIx(2)))]);
    xlabel('t'); ylabel(vNames{varIx});
    
end

end