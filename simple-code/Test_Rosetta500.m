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

param.k=50;
param.epsObj=1e-16;
param.lambda=.01;
param.mu=.005;
param.maxIter=20;
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
imagesc(abs(S1)>0);
axis square;
subplot(1,2,2);
imagesc(abs(M1)>0);
axis square;

% %%
% Z = linkage(M1,'ward');
% % [H,T,OUTPERM] = dendrogram(Z) ;
% [Cres,I]=order_of_tree(Z);
% figure(3);clf;
% subplot(1,2,1);
% imagesc(min(abs(S1),10));
% axis square;
% subplot(1,2,2);
% imagesc(min(abs(M1(I,I)),10));
% axis square;

%%
keyboard
[J]=grayorder(full(set1~=0));
figure(3);clf;
subplot(1,2,1);
imagesc(min(abs(S1(J,J)),10));
axis square;
subplot(1,2,2);
imagesc(min(abs(M1(J,J)),10));
axis square;

%%
[~, K] = sort(erdata);
figure(4); clf;
imagesc(Xn(:,K));
