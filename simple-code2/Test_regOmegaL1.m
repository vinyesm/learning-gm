
clear all
close all
clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/


% % example 1
% p = 20;
% k = 5;
% ix = 1:(3*k);
% jx = kron(1:3,ones(1,k));
% supp = sparse(ix, jx,true,p,3);
% atoms_u = sparse(ix, jx,randn(3*k,1),p,3);
% norm_atoms_u = sqrt(sum(atoms_u,1));
% atoms_u=bsxfun(@rdivide,atoms_u,norm_atoms_u);
% coeff = abs(randn(3,1));
% coeff_atoms_u = bsxfun(@times, sqrt(coeff)', atoms_u);
% M=full(coeff_atoms_u*coeff_atoms_u');
% S=triu(randn(p,p));
% S=S./max(abs(S(:)));
% S=(abs(S)>.8).*S;
% S=.5*(S+S');
% D=ones(p,1);
% emin=eigs(S-M,1,'sa');
% if emin<0 %GO
%     S=S-2*emin*diag(D);
% end
% X=inv(S-M)^.5; %X*M*X=X
% Y=eye(p);
% param.lambda=.01/p;
% param.mu=.01/p;
% %keyboard;


%example 2
p = 20;
k = 5;
tau=0.05;
invsigma=1;
h=3;

% Construction of a real-valued adjacency matrix A (columns index real-valued edges)
%       Ah are connections from the hidden variables to the rest
%       Al are connection into the graph from hidden variables
%       As are connection from and into the observed part

% As
C=sparse(double(rand(p)>(1-tau))); %sparsity pattern for sparse part
[Is,Js,vs]=find(C);
m=size(Is,1);
As=0.5*(sparse(Is,(1:m)',randn(m,1),p,m)+sparse(Js,(1:m)',randn(m,1),p,m));
% Ah & Al
ix = 1:(h*k);
jx = kron(1:h,ones(1,k));
supp = sparse(ix, jx,true,p,h);
Ah = sparse(ix, jx,sign((rand(h*k,1)-0.5)),h*k,h)';
Al=[diag(sign((rand(h*k,1)-0.5)));zeros(p-h*k,h*k)];
A=[Ah,zeros(h,m);Al,As]; 
K=A*A'; % pre-precision matrix
K=K+invsigma.^2*eye(h+p); %precision matrix for the entire graph
S=K(h+(1:p),h+(1:p)); %sparse part of the precision matrix
S=0.5*(S+S'); %stabilizing numerically
nnzS=sum(sum((S-diag(diag(S))~=0)));
fprintf('S has %d off-diagonal non-zeroes\n',nnzS);
M=K(h+(1:p),1:h)*(K(1:h,1:h)\(K(h+(1:p),1:h)')); %=K_ho (K_oo)^-1 K_oh
M=0.5*(M+M'); %stabilizing numerically
X=inv(S-M)^.5; %X*M*X=X
Y=eye(p);
param.lambda=.5;
param.mu=.005;

pause()
%keyboard

inputData.X = X;
inputData.Y = Y;
param.k=k;
set = supp;
param.epsStop=1e-4;
param.maxIter=100;
param.maxNbAtoms=1000;
param.verbose=2;


[ output, hist ] = regOmegaL1( inputData, param, set);


figure(1);clf
semilogy(hist.reldgl1,'r');hold on
semilogy(hist.reldgom,'b');
legend('rel dg l1','rel dg om');
title('relative duality gaps of the subproblems');

figure(2);clf
semilogy(hist.reldg,'k');
title('relative global duality gap');

figure(3);clf
semilogy(hist.objective,'k');
title('objective');

%dif=sign(output.S)+sign(S);
%full(output.atoms_u)
