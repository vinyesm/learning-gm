function [dualityGap] = get_dg_omega_asqp(L,Z1,param,ActiveSet,inputData,obj,loss,pen)

M=L+Z1;
H = gradient(M,inputData,param); % gradient

temp_om = -1;
for i = 1:length(ActiveSet.I)
    currS=-H(ActiveSet.I{i}, ActiveSet.I{i});
    currS=0.5*(currS+currS');
    
    if rcond(currS)<1e-14 || sum(sum(isnan(currS)))
        fprintf('get_dg_l1_omega_asqp: bad conditioned or nan\n');
    end
    
    currTemp = eigs(currS,1,'la');
%     currTemp = currTemp/(param.cardfun(ActiveSet.k{i}));
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

% duality gap

kappa = shrink * (inputData.X1*M*inputData.X2-inputData.Y);
G=inputData.X1*kappa*inputData.X2;
U=inputData.X1*L*inputData.X2;
dgf=.5*norm(U - kappa -inputData.Y,'fro')^2;
dgL= pen + trace(L*G);
fprintf('dgf=%f  dgL=%f  omega_pol=%f  (<1)\n',dgf,dgL,temp_om)
dualityGap = dgf+dgL;

keyboard;


if dualityGap<0 && abs(dualityGap)>1e-10
    fprintf('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dualityGap=abs(dualityGap);
end

end

