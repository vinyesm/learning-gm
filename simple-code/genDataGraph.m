function [X] =  genDataGraph(C,n)

p= size(C,1);
invsigma =1;

%C=sparse(double(rand(p)>(1-tau))); %sparsity pattern for sparse part
[Is,Js,vs]=find(C);
m=size(Is,1);
As=0.5*(sparse(Is,(1:m)',randn(m,1),p,m)+sparse(Js,(1:m)',randn(m,1),p,m));
% % Ah & Al
% ix = 1:(h*k);
% jx = kron(1:h,ones(1,k));
% supp = sparse(ix, jx,true,p,h);
% Ah = sparse(ix, jx,sign((rand(h*k,1)-0.5)),h*k,h)';
% Al=[diag(sign((rand(h*k,1)-0.5)));zeros(p-h*k,h*k)];
% A=[Ah,zeros(h,m);Al,As]; 
% K=A*A'; % pre-precision matrix

Ks=As*As'; % pre-precision matrix
Ks=Ks+invsigma.^2*eye(p); %precision matrix for the entire graph
S=Ks; %sparse part of the precision matrix
S=0.5*(S+S'); %stabilizing numerically
nnzS=sum(sum((S-diag(diag(S))~=0)));
fprintf('S has %d off-diagonal non-zeroes\n',nnzS);

%M=K(h+(1:p),1:h)*(K(1:h,1:h)\(K(h+(1:p),1:h)')); %=K_ho (K_oo)^-1 K_oh
%M=0.5*(M+M'); %stabilizing numerically

M=zeros(p);
X0=inv(S-M)^.5; %X*M*X=X
X=randn(n,p)*X0;
X=((1./n)*(X'*X));
