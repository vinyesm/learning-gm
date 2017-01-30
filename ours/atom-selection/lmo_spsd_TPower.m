function [uBest,kBest,allVal] = lmo_spsd_TPower(A,param)

B=0.5*(A+A');
% emin=eigs(A,1,'sa');
% if emin>0
%     emin=0;
% end
% B=A-1.1*emin*eye(size(A,1));
lambdaBest= -inf;
kBest=0;
allVal=zeros(size(A,1),1);

for k=1:size(A,1)
    if (param.cardfun(k) ~= inf)
        cf=param.cardfun(k);
        options.verbose=0;
        options.optTol=1e-8;
        options.maxIter=1000;
        options.cardinality_vec=k;
        [u,lambda] = TPower_SPCA(A, options);
        lambda=lambda/cf;
        allVal(k)=lambda;
        if lambdaBest < lambda
            uBest = u/sqrt(cf);
            lambdaBest = lambda;
            cfBest=cf;
            kBest=k;
        end
    end
end

lambdaBest=lambdaBest*cfBest;
if lambdaBest<0
    fprintf('in lmo_spca no descent direction\n');
    keyboard
end

end