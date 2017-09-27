close all; clc; clear all;
%% if using chandrasekaran method
% HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
 HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

%% example
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

%% param
Sigma = X;
param.k=k;

param.epsObj=1e-16;
param.lambda=.005;
param.mu=.005;
param.maxIter=50;
param.maxNbBlocks=100;
param.verbose=2;


%% logdetOmegaL1 initialised with true support
[S1,M1,L1,U1,hist_ch1,set1] = logdetOmegaL1(Sigma,param,supp);
K1=S1-M1;
%log(det(K1))-trace(K1*Sigma);


figure(1);clf
semilogy(hist_ch1.objective,'k');
title(['objective logdetOmegaL1 initialised with true support fend=' num2str(hist_ch1.objective(end))]);

figure(2);clf;
subplot(2,2,1);
imagesc(abs(S));
axis square;
subplot(2,2,2);
imagesc(abs(M));
axis square;
subplot(2,2,3);
imagesc(abs(S1));
axis square;
subplot(2,2,4);
imagesc(abs(M1));
axis square;

%% logdetOmegaL1 with no starting support
[S2,M2,L2,U2,hist_ch2,set2] = logdetOmegaL1(Sigma,param,inf);
K2=S2-M2;
%log(det(K2))-trace(K2*Sigma);

figure(3);clf
semilogy(hist_ch2.objective,'k');
title(['objective logdetOmegaL1 with no starting support fend=' num2str(hist_ch2.objective(end))]);

figure(4);clf;
subplot(2,2,1);
imagesc(abs(S));
axis square;
subplot(2,2,2);
imagesc(abs(M));
axis square;
subplot(2,2,3);
imagesc(abs(S1));
axis square;
subplot(2,2,4);
imagesc(abs(M1));
axis square;