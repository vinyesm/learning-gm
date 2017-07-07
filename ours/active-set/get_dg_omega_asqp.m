function [dualityGap] = get_dg_omega_asqp(L,Z1,param,ActiveSet,inputData,obj,loss,pen)

H = gradient(L+Z1,inputData,param); % gradient
Y=inputData.Y;
inputData.Y=inputData.Y+inputData.X1*(-Z1)*inputData.X2;

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

%duality gap

kappa = shrink * (inputData.X1*L*inputData.X2-inputData.Y);
% U=inputData.Y+kappa;
% dg1=.5*norm(U-inputData.Y,'fro')^2-trace(U*kappa);
% dg2= pen + trace(L*inputData.X1*kappa*inputData.X2);
U=inputData.X1*L*inputData.X2;
dg1=.5*norm(U - kappa -inputData.Y,'fro')^2;
dg2= pen + trace(L*inputData.X1*kappa*inputData.X2);
% M=inputData.X1*kappa*inputData.X2;
% M=.5*(M+M');
% ll=eigs(-M(ActiveSet.I{i}, ActiveSet.I{i})/param.lambda,1,'la');
fprintf('dg1=%f  dg2=%f  omega_pol=%f  (<1)\n',dg1,dg2,temp_om*shrink)
dualityGap = dg1+dg2;

% keyboard;

inputData.Y=Y;


if dualityGap<0 && abs(dualityGap)>1e-10
    fprintf('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dualityGap=abs(dualityGap);
end

end

