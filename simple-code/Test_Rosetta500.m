clear all;close all;clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/

HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
% HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

load('../genedata/BC.mat')


%% param

param.k=50;
param.epsObj=1e-16;
param.lambda=.005;
param.mu=.005;
param.maxIter=5;
param.maxNbBlocks=100;
param.verbose=2;


%% logdetOmegaL1 initialised with true support

[S1,M1,L1,U1,hist_ch1,set1] = logdetOmegaL1(Sigma,param,inf);
K1=S1-M1;
%log(det(K1))-trace(K1*Sigma);


figure(1);clf
semilogy(hist_ch1.objective,'k');
title(['objective logdetOmegaL1 initialised with true support fend=' num2str(hist_ch1.objective(end))]);

%%



figure(2);clf;
subplot(1,2,1);
imagesc(abs(S1));
axis square;
subplot(1,2,2);
imagesc(abs(M1)>0);
axis square;
