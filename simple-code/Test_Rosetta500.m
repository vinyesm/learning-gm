
clear all
close all
clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/

load('../genedata/BC.mat')
% example 1
p = size(Sigma,1);
k = 30;
Y=eye(p);
inputData.X = Sigma;
inputData.Y = Y;
param.k=k;
param.epsStop=1e-4;
param.lambda=.0000001/p*1000;
param.mu=.001/p*1000;
param.maxIter=250;
param.maxNbAtoms=1000;
param.verbose=2;

%%
[ output, hist ] = regOmegaL1( inputData, param, inf);


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

%%
S0 = output.S;
L0 = output.M;
figure(4);clf
subplot(1,2,1)
imagesc(abs(S0-diag(diag(S0)))); colormap hot;
axis square
subplot(1,2,2)
imagesc(abs(L0)>0); colormap hot;
axis square

%dif=sign(output.S)+sign(S);
%full(output.atoms_u)
