function [uBest,kBest,lambdaBest] = lmo_spsd_TPower(A,param)
both=0;
B=0.5*(A+A');

% emin=eigs(B,1,'sa');
emin=-norm(B,'fro')/2;

if emin<0
    B=B-2*emin*eye(size(A,1));    
end
% emin=eigs(A,1,'sa');
% if emin>0
%     emin=0;
% end
% B=A-1.1*emin*eye(size(A,1));
lambdaBest= -inf;
lambda_am = -inf;
kBest=0;
allVal=zeros(size(B,1),1);

for k=1:size(B,1)
    if (param.cardfun(k) ~= inf)
        cf=param.cardfun(k);
        options.verbose=0;
        options.optTol=1e-8;
        options.maxIter=1000;
        options.cardinality_vec=k;
        options.initType=1; %default 2
        [uTP,lambdaTP] = TPower_SPCA(B, options);
        if both
            [u ,nIter ,~, ~] = spca_am(rand(size(A,1),1), B, 2, 0, 0,k, 1000, 1e-6);
            lambda_am=u'*B*u;
            lambda=lambda_am;
            fprintf('TPower : %f am : %f\n',lambdaTP,lambda_am);
        end
        if lambdaTP>lambda_am
            lambda=lambdaTP;
            u=uTP;
        end
        if emin<0
            lambda=(lambda+2*emin)/cf;
        else
            lambda=lambda/cf;
        end
        allVal(k)=lambda;
        if  lambdaBest < lambda
            uBest = u/sqrt(cf);
            lambdaBest = lambda;
            cfBest=cf;
            kBest=k;
        end
    end
end
% if emin<0
%     lambdaBest=lambdaBest*cfBest+2*emin;   
% end
% lambdaBest=lambdaBest/cfBest;

if lambdaBest<0
    fprintf('in lmo_spca no descent direction\n');
%     keyboard;
end

end