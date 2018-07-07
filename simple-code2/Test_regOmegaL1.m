
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
% param.lambda=.01/p*1000;
% param.mu=.01/p*1000;
% %keyboard;


%example 2
p = 160;
n= 2000;
k = 35;
tau=0.01;
invsigma=1;
h=4;

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
X0=inv(S-M)^.5; %X*M*X=X
X=randn(n,p)*X0;
X=((1./n)*(X'*X))^.5;
Y=eye(p);
% param.lambda=.02;
% param.mu=.002;
param.lambda=.02;
param.mu=.019;
param.epsStop=1e-4;

%pause()
%keyboard

inputData.X = X;
inputData.Y = Y;
param.k=k;
set = supp;
%param.epsStop=1e-3;
param.maxIter=100;
param.maxNbAtoms=1000;
param.verbose=2;


[ output, hist ] = regOmegaL1( inputData, param, set );


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

figure(17)
imagesc(output.M); colormap gray; axis square
figure(18)
imagesc(M);axis square

figure(19)
imagesc(abs(-full(output.S))>1e-6); colormap gray; axis square
figure(20)
imagesc(S~=0);axis square

if 0,
figure(10);clf
imagesc(-abs(((X'*X)))); colormap gray; axis square
end
%%
greenColorMap = [ones(1, 132), linspace(1, 0, 124)];
redColorMap = ones(1,256);%[linspace(0, 1, 124), ones(1, 132)];
blueColorMap = [linspace(0, 1, 62), ones(1,132), linspace(1, 0, 62)];
colorMap = [redColorMap; greenColorMap; blueColorMap]';
% Apply the colormap.
colormap(colorMap);
%colormap(flipud(gray))
if 0
mkdir('fig');
figure(21)
imagesc(output.M); colormap(colorMap); axis square
%set(gca,'FontSize',40)
colorbar('southoutside')
print('fig/outputM','-depsc')

figure(22)
%set(gca,'FontSize',18)
imagesc(M); colormap(colorMap); axis square
colorbar('southoutside')
print('fig/M','-depsc')
end
