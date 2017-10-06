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
%run('../DREAM5/Dream/pp_net1');
load('../DREAM5/Dream/pp_net1.mat')

%% param

param.k=15;
param.epsObj=1e-16;
param.lambda=.004;
param.mu=.002;
param.maxIter=10;
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

[J]=grayorder(full(set1~=0));
figure(3);clf;
subplot(2,2,1);
imagesc(min(abs(S1),10));
axis square;
subplot(2,2,2);
imagesc(min(abs(M1),10));
axis square;
subplot(2,2,3);
imagesc(min(abs(S1(J,J)),10));
axis square;
subplot(2,2,4);
imagesc(min(abs(M1(J,J)),10));
axis square;

figure(4);clf;
subplot(2,3,1);
imagesc(abs(S1)>1e-2); colormap gray
axis square;
subplot(2,3,2);
imagesc(abs(M1)>1e-2);
axis square;
subplot(2,3,3);
imagesc(Ac);
axis square;
subplot(2,3,4);
imagesc(abs(S1(J,J))>1e-2);
axis square;
subplot(2,3,5);
imagesc(abs(M1(J,J))>1e-2);
axis square;
subplot(2,3,6);
imagesc(Ac(J,J));
axis square;



%%
% [~, K] = sort(erdata);
% figure(4); clf;
% imagesc(Xn(:,K));
