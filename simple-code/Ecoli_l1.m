clear all;close all;clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/

%HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
 HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

addpath ../reorder/

%load('../genedata/BC.mat')
run('../DREAM5/Dream/pp_net1');

%% param

param.k=30;
param.epsObj=1e-16;
param.lambda=1000;
param.maxIter=50;
param.maxNbBlocks=100;
param.verbose=2;

pos = sum(sum((Ac==1)));
neg = sum(sum((Ac==0)));
tpr = [];
fpr = [];

%% logdetOmegaL1 initialised with true support
mus = 2.^linspace(-16,1,30);
for mu=mus
    param.mu=mu;
    [S1,M1,L1,U1,hist_ch1,set1] = logdetOmegaL1(Sigma,param,inf);
    K1=S1-M1;

    Sc = abs(S1-diag(diag(S1)))>0;
    tp = sum(sum((Ac==1) & Sc));
    fp = sum(sum((Ac==0) & Sc));
    tpr = [tpr, tp/pos];
    fpr = [fpr, fp/(fp+neg)];
end


figure(1);clf
plot(fpr,tpr,'.');hold on;
plot([0 1],[0,1],'r-'); 
axis([0 1 0 1])
