clear all
close all
clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/

% % example 1
% p = 6;
% k = 2;
% ix = 1:(3*k);
% jx = kron(1:3,ones(1,k));
% supp = sparse(ix, jx,true,p,3);
% atoms_u = sparse(ix, jx,randn(3*k,1),p,3);
% atoms_u = normc(atoms_u);
% coeff = abs(randn(3,1));
% coeff_atoms_u = bsxfun(@times, sqrt(coeff)', atoms_u);
% M=full(coeff_atoms_u*coeff_atoms_u');
% X=pinv(M); %X*M*X=X
% Y=X;
% 
% inputData.X = X;
% inputData.Y = Y;
% param.k=k;
% set = supp;
% param.eps=1e-6;
% param.lambda=.1;
% param.maxIter=50;
% param.maxNbAtoms=20;
% 
% [ output1, hist1 ] = regOmegaRestricted( inputData, param, set);
% keyboard;

%%
% example 2
p = 1000;
k = 10;
ix = [1:10 6:15 11:20];
jx = kron(1:3,ones(1,10));
supp = sparse(ix, jx,true,p,3);
%atoms_u = sparse(ix, jx,randn(3*k,1),p,3);
atoms_u = sparse(ix, jx,ones(3*k,1),p,3);
atoms_u = normc(atoms_u);
coeff = abs(randn(3,1));

coeff_atoms_u = bsxfun(@times, sqrt(coeff)', atoms_u);
M=full(coeff_atoms_u*coeff_atoms_u');
X=pinv(M); %X*M*X=X
Y=X;

inputData.X = X;
inputData.Y = Y;
param.k=k;
set = supp;
param.eps=1e-6;
param.lambda=1/p;
param.maxIter=50;
param.maxNbAtoms=20;
param.verbose=2;

ntest=1;
objTest=zeros(1,ntest);
for i=1:ntest
[ output2, hist2 ] = regOmegaRestricted( inputData, param, set);
objTest(i)=hist2.objective(end);
end
%keyboard;
outa = output2.atoms_u(:,1:output2.atomCount);
% full(output2.atoms_u(:,1:output2.atomCount))
% full(atoms_u)

figure(1)
semilogy(hist2.dualityGap);

figure(2);clf
semilogy(hist2.objective,'k');
title('objective')

