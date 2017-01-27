function [ A,M,S,E, as,out] = sparse_omega_lgm( inputData, param)
% min_(A,M,S) .5|C^.5*A*C^.5|^2 + mu|S|_1 + lambda Omega_psd,k(M)
% s.t. A>=0 and M>=0
% using ADMM

debug=1;

options.cardfun=param.cardfun;
lossopt='quad';%options.loss;

C05=inputData.X1;
C=C05*C05;
p=size(C,1);
rho=param.rho;
%1/norm(C^4,'fro');
lambda=param.lambda;
mu=param.mu;
max_iter=200;
%init A,M,S,E
A=zeros(p);
M=zeros(p);
S=zeros(p);
E=zeros(p);
Mnorm=0;


if debug
    obj=zeros(1,max_iter*4);
    auglag=zeros(1,max_iter*4);
    residual=zeros(1,max_iter*4);
    count=0;
    obj2=zeros(1,max_iter);
    auglag2=zeros(1,max_iter);
    residual2=zeros(1,max_iter);
end

if strcmp(lossopt,'quad')
    loss=@(A) .5*norm(C05*A*C05-eye(p),'fro')^2;
elseif strcmp(lossopt,'ml')
    error('still not implemented')
end
aug_lag=@(A,M,Mnorm,S,E) loss(A)+mu*sum(abs(S(:)))+lambda*Mnorm-trace(E'*(S-M-A))+rho/2*norm(S-M-A,'fro')^2;
objective=@(A,M,Mnorm,S,E) loss(A)+mu*sum(abs(S(:)))+lambda*Mnorm;


for i=1:max_iter
   
    % A with ls eta=norm(G,'fro')^2/(norm(C^.5*G*C^.5,'fro')^2 + rho*norm(G,'fro')^2)
    G=C*A*C-C+E-rho*(S-M-A);
    eta=norm(G,'fro')^2/(norm(C^.5*G*C^.5,'fro')^2 + rho*norm(G,'fro')^2);
    A=A-eta*G;
    A=proj_psd(A);
    if debug
        count=count+1;
        obj(count)=objective(A,M,Mnorm,S,E);
        auglag(count)=aug_lag(A,M,Mnorm,S,E);
        residual(count)=norm(A+M-S,'fro')^2;
    end
    S=soft_threshold(M+A+E/rho,mu/rho);
    if debug
        count=count+1;
        obj(count)=objective(A,M,Mnorm,S,E);
        auglag(count)=aug_lag(A,M,Mnorm,S,E);
        residual(count)=norm(A+M-S,'fro')^2;
    end
    [M Mnorm as]=prox_omega(S-A-E/rho,lambda/rho,options);
    if Mnorm<0
        error('negative Mnorm\n');
    end
    if debug
        count=count+1;
        obj(count)=.5*norm(C05*A*C05-eye(p),'fro')^2+mu*sum(abs(S(:)))+lambda*Mnorm;
        auglag(count)=aug_lag(A,M,Mnorm,S,E);
        residual(count)=norm(A+M-S,'fro')^2;
    end
    S=soft_threshold(M+A+E/rho,mu/rho);
    if debug
        count=count+1;
        obj(count)=objective(A,M,Mnorm,S,E);
        auglag(count)=aug_lag(A,M,Mnorm,S,E);
        residual(count)=norm(A+M-S,'fro')^2;
    end
    E=E-rho*(S-M-A);
    if debug
        count=count+1;
        obj(count)=objective(A,M,Mnorm,S,E);
        auglag(count)=aug_lag(A,M,Mnorm,S,E);
        residual(count)=norm(A+M-S,'fro')^2;
    end
    
    
%     if debug
%         fprintf(['------------------------------------------------------------\n']);
%         fprintf(['sparse+omega iter    ' num2str(i) '\n']);
%         fprintf(['------------------------------------------------------------\n']);
%         obj(i)=.5*norm(C05*A*C05-eye(p),'fro')^2+mu*sum(abs(S(:)))+lambda*Mnorm;
%         residual(i)=norm(A+M-S,'fro')^2;
%     end
    if debug
        obj2(i)=objective(A,M,Mnorm,S,E);
        auglag2(i)=aug_lag(A,M,Mnorm,S,E);
        residual2(i)=norm(A+M-S,'fro')^2;
        out.obj=objective(S-M,M,Mnorm,S,E);
        out.Mnorm=Mnorm;
    end
    
end

if debug
    figure(20);clf;
    subplot(2,3,1);
    semilogy(obj,'.');
    title('objective');
    subplot(2,3,2);
    semilogy(auglag,'.');
    title('auglag');
    subplot(2,3,3);
    semilogy(residual,'.');
    title('residual');
    subplot(2,3,4);
    semilogy(obj2,'r.');
    title('objective');
    subplot(2,3,5);
    semilogy(auglag2,'r.');
    title('auglag');
    subplot(2,3,6);
    semilogy(residual2,'r.');
    title('residual');
    keyboard;
end

end

