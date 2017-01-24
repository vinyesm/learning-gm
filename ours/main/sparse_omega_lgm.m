function [ A,M,S,E] = sparse_omega_lgm( inputData, param)
% min_(A,M,S) .5|C^.5*A*C^.5|^2 + mu|S|_1 + lambda Omega_psd,k(M)
% s.t. A>=0 and M>=0
% using ADMM

debug=1;

options.cardfun=param.cardfun;

C05=inputData.X1;
C=C05*C05;
p=size(C,1);
rho=param.rho;
%1/norm(C^4,'fro');
lambda=param.lambda;
mu=param.mu;
max_iter=500;
%init A,M,S,E
A=zeros(p);
M=zeros(p);
S=zeros(p);
E=zeros(p);
Mnorm=0;


if debug
    obj=zeros(1,max_iter);
    eq_constraint=zeros(1,max_iter);
end

aug_lag=@(A,M,Mnorm,S,E) .5*norm(C05*A*C05-eye(p),'fro')^2+mu*sum(abs(S(:)))+lambda*Mnorm-trace(E'*(S-M-A))+rho/2*norm(S-M-A,'fro')^2;

for i=1:max_iter
% A with ls eta=norm(G,'fro')^2/(norm(C^.5*G*C^.5,'fro')^2 + rho*norm(G,'fro')^2)
fprintf(['0) augmented lagrangian  ' num2str(aug_lag(A,M,Mnorm,S,E)) '\n']);
G=C*A*C-C+E-rho*(S-M-A);
eta=norm(G,'fro')^2/(norm(C^.5*G*C^.5,'fro')^2 + rho*norm(G,'fro')^2);
A=A-eta*G;
A=proj_psd(A);
fprintf(['1) augmented lagrangian  ' num2str(aug_lag(A,M,Mnorm,S,E)) '\n']);
S=soft_threshold(M+A+E/rho,mu/rho);
fprintf(['2) augmented lagrangian  ' num2str(aug_lag(A,M,Mnorm,S,E)) '\n']);
[M Mnorm]=proj_omega(S-A-E/rho,lambda/rho,options);
fprintf(['3) augmented lagrangian  ' num2str(aug_lag(A,M,Mnorm,S,E)) '\n']);
E=E-rho*(S-M-A);
fprintf(['4) augmented lagrangian  ' num2str(aug_lag(A,M,Mnorm,S,E)) '\n']);

if debug
    fprintf(['------------------------------------------------------------\n']);
    fprintf(['sparse+omega iter    ' num2str(i) '\n']);
    fprintf(['------------------------------------------------------------\n']);
    obj(i)=.5*norm(C05*A*C05-eye(p),'fro')^2+mu*sum(abs(S(:)))+lambda*Mnorm;
    eq_constraint(i)=norm(A+M-S,'fro')^2;
end

end

if debug
    figure(20);clf;
    subplot(1,2,1);
    plot(obj,'.');
    title('objective');
    subplot(1,2,2)
    plot(eq_constraint,'.');
    title('eq constraint');  
    keyboard;
end

end

