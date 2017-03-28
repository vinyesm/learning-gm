function [dualityGap] = get_dg_l1_omega_asqp(Z,param,ActiveSet,inputData,obj,loss,pen)

H = gradient(Z,inputData,param); % gradient
temp_l1=max(abs(H(:)));

temp_om = -1;
if param.PSD
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
else
    for i = 1:length(ActiveSet.I)
        currTemp = norm(H(ActiveSet.I{i}, ActiveSet.J{i}));
        if currTemp>temp_om
            temp_om = currTemp;
        end
    end
end

%test if temp_om<0 constraint already satisfied. We want to select 1 in min(1, param.lambda / temp_om)
if temp_om<0
    shrink=min(1,param.mu / temp_l1);
else
    shrink=min(1, min(param.lambda / temp_om, param.mu / temp_l1));
end

%duality gap

kappa = shrink * (inputData.X1*Z*inputData.X2-inputData.Y);
dualityGap = obj + sum(sum( inputData.Y .* kappa )) + 1 / 2 * sum(sum( kappa .* kappa ));
gapLoss = loss + sum(sum( inputData.Y .* kappa )) + 1 / 2 * sum(sum( kappa .* kappa )) -sum(sum((inputData.X1*Z*inputData.X2).*kappa));
gapPen = pen + sum(sum((inputData.X1*Z*inputData.X2).*kappa));


if dualityGap<0 && abs(dualityGap)>1e-10
    fprintf('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dualityGap=abs(dualityGap);
end

end

