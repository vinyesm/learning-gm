function [dg1, dg2]= compute_dg_S_L(S,L,inputData,param,ActiveSet)

%computing dg
M=L+S;
H = inputData.X1'*(inputData.X1*(S+L)*inputData.X2 - inputData.Y)*inputData.X2';


%% dg1 (S)
temp_l1 = max(abs(H(:)));
shrink=min(1,param.mu / temp_l1);
kappa = shrink * (inputData.X1*M*inputData.X2-inputData.Y);
G=inputData.X1*kappa*inputData.X2;
dg_f=.5*norm(inputData.X1*M*inputData.X2 - kappa -inputData.Y,'fro')^2;
dg_S=param.mu*sum(abs(S(:)))+trace(G*S);
dg1=dg_f+dg_S;

fprintf('    dgf=%f  dgS=%f  dg1=%f\n',dg_f,dg_S,dg1)

if dg_f<0 || dg_S<0  && abs(dg1)>1e-10
    fprintf('Negative duality gap\n');
    keyboard;
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dg1=abs(dg1);
end

%% dg2 (L)
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

if temp_om<0
    shrink=1;
else
    shrink=min(1, param.lambda / temp_om);
end

% duality gap

kappa = shrink * (inputData.X1*M*inputData.X2-inputData.Y);
G=inputData.X1*kappa*inputData.X2;
U=inputData.X1*M*inputData.X2;
dg_f=.5*norm(U - kappa -inputData.Y,'fro')^2;
dg_L= param.lambda*sum(ActiveSet.alpha(1:ActiveSet.atom_count)) + trace(L*G);
dg2 = dg_f+dg_L;
fprintf('    dgf=%f  dgL=%f  dg2=%f\n',dg_f,dg_L,dg2)

if dg_f<0 || dg_L<0 && abs(dg2)>1e-10
    fprintf('Negative duality gap\n');
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dg2=abs(dg2);
end
