function [uBest,kBest] = lmo_spca(A,param)
B=0.5*(A+A');
% emin=eigs(A,1,'sa');
% if emin>0
%     emin=0;
% end
% B=A-1.1*emin*eye(size(A,1));
lambdaBest= -inf;
kBest=0;
for k=1:size(A,1)
    if (param.cardfun(k) ~= inf)
        cf=param.cardfun(k);
        param.k=k;
        param.q=k;
        [u,lambda] = TPI(B,param);
        lambda=(u'*A*u)/cf;
        if lambdaBest < lambda
            uBest = u/sqrt(cf);
            lambdaBest = lambda;
            kBest=k;
        end
    end
end

if lambdaBest<0
    fprintf('in lmo_spca no descent direction\n');
end

end