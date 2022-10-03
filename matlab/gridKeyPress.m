function gridKeyPress(src, evt, clo, hi, Grid, feature)

plb=[clo.prob.par.lb];
pub=[clo.prob.par.ub];
initBounds = [plb(Grid.xix),pub(Grid.xix),plb(Grid.yix),pub(Grid.yix)];

disp(evt.Key)
switch(evt.Key)
    case 'c' %'continue' - features without initialization
        tic
        clo.shiftX0();
        clo.features();
        toc
        F=clo.getF();
%         F(F<-1e10|F>1e10)=nan;
        hi.CData=reshape(feature.fun(F),Grid.dim);
        hcb=colorbar('northoutside');
        title(hcb,feature.name,'Interpreter','none')
%         axis tight
        
    case 'g' %'go' - features with initialization
        tic
        clo.shiftX0();
        clo.features(1);
        toc
        F=clo.getF();
%         F(F<-1e10|F>1e10)=nan;
        hi.CData=reshape(feature.fun(F),Grid.dim);
        hcb=colorbar('northoutside');
        title(hcb,feature.name,'Interpreter','none')
%         axis tight
        
    case 't' %'transient'
        tic
        clo.shiftX0();
        clo.transient();
        toc
        Xf=clo.getXf();
        hi.CData=reshape(Xf(:,clo.op.fVarIx),Grid.dim);
        hcb=colorbar('northoutside');
        title(hcb,clo.prob.varNames(clo.op.fVarIx),'Interpreter','none')
%         axis tight
        
    case 'f' %select feature to plot %%%TODO: add computed features
%         featureSelectDialog(feature.fun,feature.name)
%         hi.CData=reshape(feature.fun(F),grid.dim);
%         hcb=colorbar('northoutside');
%         title(hcb,feature.name)
        
    case 'p' %change parameter values
        %select p1, p2?
        %new plb1 plb2 pub1 pub2
        %generate new P(:)
        %upload to GPU
        %run, or wait for 'i'/'c'?
        
    case 'i' %set initial condition data
        
    case 'r' %randomize initial condition
        
        x0lb=[clo.prob.var.lb];
        x0ub=[clo.prob.var.ub];
        X0=x0lb+rand(clo.nPts,length(x0lb)).*(x0ub-x0lb);
        clo.setX0(X0);
        
%         tic
%         clo.transient();
%         toc
%         Xf=clo.getXf();
        hi.CData=reshape(X0(:,clo.op.fVarIx),Grid.dim);
        hcb=colorbar('northoutside');
        title(hcb,clo.prob.varNames(clo.op.fVarIx),'Interpreter','none')
%         axis tight
        
    case 'z' %zoom using rbbox
        
        tmp=get(hi.Parent,'ButtonDownFcn'); %
        set(hi.Parent,'ButtonDownFcn',[]); %temporarily turn off clicktrajectory
        
        gridAxis=gca;
        
        k=waitforbuttonpress;
        if k==0
            point1 = gridAxis.CurrentPoint;    % location when mouse clicked
            rbbox;                 % rubberbox
            point2 = gridAxis.CurrentPoint;    % location when mouse released
            point1 = point1(1,1:2);            % extract x and y
            point2 = point2(1,1:2);
            lb = min(point1,point2);           % calculate bounds
            ub = max(point1,point2);           % calculate bounds

            if all(ub-lb>0)
                setBounds(lb,ub);
            else
                disp('Try selecting new bounds more slowly...')
            end
        else
            disp('Window Zoom-in canceled by keypress')
        end
        
        set(hi.Parent,'ButtonDownFcn',tmp); %restore clicktrajectory
    
    case 'add'
        
    case 'subtract'
        
    case '1'
        setBounds(initBounds([1,3]),initBounds([2,4]));
        
end  

    function setBounds(lb,ub)
        p1=linspace(lb(1),ub(1),Grid.dim(1));
        p2=linspace(lb(2),ub(2),Grid.dim(2));
        [P1,P2]=meshgrid(p1,p2);
        P=clo.P;
        P(:,Grid.xix)=P1(:);
        P(:,Grid.yix)=P2(:);
        clo.setP(P);
        axis(gca,[lb(1),ub(1),lb(2),ub(2)]);
        hi.XData=[lb(1),ub(1)];
        hi.YData=[lb(2),ub(2)];
        hi.CData=nan(size(hi.CData));
    end

    function featureSelectDialog(feature)
    %     %changes the feature.fun and feature.name
    %     hfs=figure('Name','Select Feature to Plot','WindowStyle','modal');
    %     hfl=uicontrol('Style','popupmenu','String',clo.fNames);
    %     hfb=uibutton('Text','OK','ButtonPushedFcn',selectFeature);
    % 
    %     
    %     function selectFeature(src,evt)
    %         feature.fun=@
    %     end
    end


end

function parameterDialog()
end
