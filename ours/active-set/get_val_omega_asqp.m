function [loss,pen,obj,dualityGap,time]=get_val_omega_asqp(Z,Z1,L,ActiveSet,inputData,param)

cf=zeros(ActiveSet.atom_count,min(ActiveSet.atom_count,1));
card=sum(ActiveSet.atoms(:,1:ActiveSet.atom_count)~=0);
for i=1:ActiveSet.atom_count
    cf(i)=min(param.cardfun(card(i):end));
end

if param.f==4
    currloss = .5*norm(inputData.Y - inputData.X1*Z*inputData.X2, 'fro')^2;
    loss= currloss;
    if ~isempty(ActiveSet.alpha)
        pen = param.lambda*dot(cf,ActiveSet.alpha);
    else
        pen = 0;
    end
    obj = currloss + pen ;
    dualityGap = get_dg_omega_asqp(L,Z1,param,ActiveSet,inputData,obj,loss,pen);
    dualityGap= real(dualityGap);
    time=toc;
elseif param.f==5
    currloss = .5*trace(Z^2*inputData.X)-trace(Z*inputData.Y);
    loss= currloss;
    if ~isempty(ActiveSet.alpha)
        pen = param.lambda*dot(cf,ActiveSet.alpha);
    else
        pen = 0;
    end
    obj = currloss + pen ;
    dualityGap = inf;
    time=toc;
end

end

