
clear all
close all
clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/


% example 1
p = 20;
k = 5;
ix = 1:(3*k);
jx = kron(1:3,ones(1,k));
supp = sparse(ix, jx,true,p,3);
atoms_u = sparse(ix, jx,randn(3*k,1),p,3);
atoms_u = normc(atoms_u);
coeff = 5*abs(randn(3,1));
coeff_atoms_u = bsxfun(@times, sqrt(coeff)', atoms_u);
M=full(coeff_atoms_u*coeff_atoms_u');
S=triu(randn(p,p));
S=S./max(abs(S(:)));
S=(abs(S)>.9).*S;
S=.5*(S+S');
D=ones(p,1);
emin=eigs(S-M,1,'sa');
S=S-2*emin*diag(D);
X=inv(S-M)^.5; %X*M*X=X
Y=eye(p);

inputData.X = X;
inputData.Y = Y;
param.k=k;
set = supp;
param.epsStop=1e-4;
param.lambda=.01/p;
param.mu=.01/p;
param.maxIter=250;
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
