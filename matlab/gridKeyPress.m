function gridKeyPress(src,evt, clo, hi, Ffun, ftitle, nGrid)

switch(evt.Key)
    case 'c' %continue
        tic
        clo.features();
        toc
        F=clo.getF();
        hi.CData=reshape(Ffun(F),nGrid);
        hcb=colorbar('northoutside');
        title(hcb,ftitle)
    case 'f' %init features
        tic
        clo.features(1);
        toc
        F=clo.getF();
        hi.CData=reshape(Ffun(F),nGrid);
        hcb=colorbar('northoutside');
        title(hcb,ftitle)
    case 't' %transient
        tic
        clo.transient();
        toc
        X0=clo.getX0();
        hi.CData=reshape(X0(:,clo.op.fVarIx),nGrid);
        hcb=colorbar('northoutside');
        title(hcb,clo.prob.varNames(clo.op.fVarIx))
    case 'p' %change parameter values
        %select p1, p2?
        %new plb1 plb2 pub1 pub2
        %generate new P(:)
        %upload to GPU
        %run, or wait for 'i'/'c'?
    case 'r' %randomize initial condition
end 

end


function [p,plu,pub]=parameterDialog()
end