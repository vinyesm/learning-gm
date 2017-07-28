function [dualityGap] = get_dg_omega_asqp(L,Z1,param,ActiveSet,inputData,obj,loss,pen)
debug=1;
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
U=inputData.X1*M*inputData.X2;
dgf=.5*norm(U - kappa -inputData.Y,'fro')^2;
dgL= pen + trace(L*G);
dualityGap = dgf+dgL;
fprintf('dgf=%f  dgL=%f  dgPS=%f  (<1)\n',dgf,dgL,dualityGap);
%%
if debug
    nb=length(ActiveSet.I);
    sumtrace=0;
    for i=1:nb
        supp=ActiveSet.I{i};
        len=length(supp);
        block=zeros(len);
        K=false(1,ActiveSet.atom_count);
        for at=1:ActiveSet.atom_count
            u=ActiveSet.atoms(:,at);
            %         keyboard;
            if len==sum(u~=0) && all(supp==find(u~=0))
                K(at)=1;
                block=block+ActiveSet.alpha(at)*(u(supp)*u(supp)');
            end
        end
        ActiveSet.block{i}=block;
        sumtrace=sumtrace+trace(ActiveSet.block{i});
        ActiveSet.supp{i}=supp;
    end
    sumtrace=param.lambda*sumtrace;
    if (sumtrace-pen)^2>1e-8
        keyboard;
    end
end
%%

if dgf<0 || dgL<0 && abs(dualityGap)>1e-10
    fprintf('Negative duality gap\n');
    %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
    dualityGap=abs(dualityGap);
end

end

