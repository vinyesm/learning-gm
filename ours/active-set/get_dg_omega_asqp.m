function [dualityGap] = get_dg_omega_asqp(Z,param,ActiveSet,inputData,obj,loss,pen)

H = gradient(Z,inputData,param); % gradient
temp_l1=max(abs(H(:)));

temp_om = -1;
for i = 1:length(ActiveSet.I)
    currS=-H(ActiveSet.I{i}, ActiveSet.I{i});
    currS=0.5*(currS+currS');
    
    if rcond(currS)<1e-14 || sum(sum(isnan(currS)))
        fprintf('get_dg_l1_omega_asqp: bad conditioned or nan\n');
    end
    
    currTemp = eigs(currS,1,'la');
    currTemp = currTemp/(param.cardfun(ActiveSet.k{i}));
    if currTemp>temp_om
        temp_om = currTemp;
    end
end

%test if temp_om<0 constraint already satisfied. We want to select 1 in min(1, param.lambda / temp_om)
if temp_om<0
    shrink=1;
else
    shrink=min(1, param.lambda / temp_om);
end

%duality gap

kappa = shrink * (inputData.X1*Z*inputData.X2-inputData.Y);
U=inputData.Y+kappa;
dg1=.5*norm(inputData.X1*U*inputData.X2-inputData.Y,'fro')^2;
dg2= pen + trace(Z,inputData.X1*kappa*inputData.X2);
dualityGap = dg1+dg2;
keyboard;


if dualityGap<0 && abs(dualityGap)>1e-10
    fprintf('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dualityGap=abs(dualityGap);
end

end

