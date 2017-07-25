function [dg_global] = get_dg_global(L,S,D,inputData,param,ActiveSet,pen)

% global duality gap after exact minimization in D
M=L+S+D;
H = gradient(M,inputData,param); % gradient
Y=inputData.Y;

[u, kBest,temp_om] = lmo_spsd_TPower(-H,param);
% temp_om = -1;
% for i = 1:length(ActiveSet.I)
%     currS=-H(ActiveSet.I{i}, ActiveSet.I{i});
%     currS=0.5*(currS+currS');
%     
%     if rcond(currS)<1e-14 || sum(sum(isnan(currS)))
%         fprintf('get_dg_l1_omega_asqp: bad conditioned or nan\n');
%     end
%     
%     currTemp = eigs(currS,1,'la');
% %     currTemp = currTemp/(param.cardfun(ActiveSet.k{i}));
%     currTemp = currTemp/(param.cardfun(ActiveSet.k{i}));
%     if currTemp>temp_om
%         temp_om = currTemp;
%     end
% end

temp_l1 = max(max(abs(H-diag(diag(H)))));

%test if temp_om<0 constraint already satisfied. We want to select 1 in min(1, param.lambda / temp_om)
shrink=min(1, param.lambda / temp_om, param.mu / temp_l1);
if shrink<0
    shrink=1;
end

%duality gap

kappa = shrink * (inputData.X1*Z*inputData.X2-inputData.Y);
G=inputData.X1*kappa*inputData.X2;
dg_f=.5*norm(inputData.X1*Z*inputData.X2 - kappa -inputData.Y,'fro')^2;
dg_S=param.mu*sum(abs(S(:)))+trace(G*S);
dg_L=param.lambda*sum(ActiveSet.alpha(1:ActiveSet.atom_count))+trace(G*L);
dg_D=trace(D*G);
dg_global=dg_f+dg_S+dg_L+dg_D;
fprintf('dg_f=%f  dg_S=%f  dg_L=%f  dg_D=%f dg_global=%f\n',dg_f,dg_S,dg_L,dg_D,dg_global);
dualityGap = dg1+dg2;

% keyboard;

inputData.Y=Y;


if dualityGap<0 && abs(dualityGap)>1e-10
    fprintf('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dualityGap=abs(dualityGap);
end

end

