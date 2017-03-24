function [loss,pen,obj,dualityGap,time]=get_val_l1_omega_asqp(Z,ActiveSet,inputData,param,cardVal)

if param.f==4
    currloss = .5*norm(inputData.Y - inputData.X1*Z*inputData.X2, 'fro')^2;
    loss= currloss;
    pen= param.lambda*sum(ActiveSet.alpha)+param.mu*sum(ActiveSet.beta);
    obj = currloss + pen ;
    dualityGap = get_dg_l1_omega_asqp(Z,param,ActiveSet,inputData,obj,loss,pen);
    dualityGap= real(dualityGap);
    time=toc;
elseif param.f==5
    currloss = .5*trace(Z^2*inputData.X)-trace(Z*inputData.Y);
    loss= currloss;
    pen= param.lambda*sum(ActiveSet.alpha)+param.mu*sum(ActiveSet.beta);
    obj = currloss + pen ;
    dualityGap = inf;
    time=toc;
end

end

