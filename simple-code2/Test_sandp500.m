clc; clear all; close all;
%load('../DREAM5/Dream/pp_net1.mat')
load('../sandp500/sandp500.mat')
addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/
addpath ../reorder

param.k=150;
param.epsObj=1e-16;
% param.lambda=.03;
% param.mu=.02;
param.lambda=.5;
param.mu=.01;
param.maxIter=5000;
param.maxNbAtoms=1000;
param.verbose=1;
param.epsStop=1e-3;

p = size(Sigma,1);
inputData.X = (Sigma+eye(p)*1e-5)^0.5;
inputData.X = 0.5*(inputData.X+inputData.X');
inputData.Y = eye(p);

% keyboard;

[ output1, hist1 ] = regOmegaL1( inputData, param, inf );

%%

figure(1);clf
imagesc(abs(Sigma));

figure(10);clf
semilogy(hist1.reldgl1,'r');hold on
semilogy(hist1.reldgom,'b');
legend('rel dg l1','rel dg om');
title('relative duality gaps of the subproblems');

figure(11);clf
semilogy(hist1.reldg,'k');
title('relative global duality gap');

figure(12);clf
semilogy(hist1.objective,'k');
title('objective');

figure(17)
imagesc(output1.M);
% figure(18)
% imagesc(M);

figure(19)
imagesc(-full(output1.S));

tt = 1;



figure(20)
subplot(3,2,1)
imagesc(min(abs(output1.S),tt));axis square; colormap jet
subplot(3,2,2)
imagesc(min(abs(output1.M),tt));axis square
subplot(3,2,3)
imagesc(min(abs(output1.S(I,I)),tt));axis square
subplot(3,2,4)
imagesc(min(abs(output1.M(I,I)),tt));axis square
subplot(3,2,5)
imagesc(min(abs(output1.S(K,K)),tt));axis square
subplot(3,2,6)
imagesc(min(abs(output1.M(K,K)),tt));axis square


%%
ua = unique(abs(output1.atoms_u')>1e-10, 'rows');
set = ua';
figure(40); imagesc(set); colormap gray

%%
[J]=grayorder(full(set~=0));
figure(21);clf;
subplot(2,2,1);
imagesc(min(abs(output1.S(K,K)),10));
axis square;
subplot(2,2,2);
imagesc(min(abs(output1.M(K,K)),10));
axis square;
subplot(2,2,3);
imagesc(min(abs(output1.S(J,J)),10));
axis square;
subplot(2,2,4);
imagesc(min(abs(output1.M(J,J)),10));
axis square;

figure(22)
subplot(3,2,1)
imagesc(min(abs(output1.S),tt)>0);axis square; colormap gray
subplot(3,2,2)
imagesc(min(abs(output1.M),tt)>0);axis square
subplot(3,2,3)
imagesc(min(abs(output1.S(I,I)),tt)>0);axis square
subplot(3,2,4)
imagesc(min(abs(output1.M(I,I)),tt)>0);axis square
subplot(3,2,5)
imagesc(min(abs(output1.S(K,K)),tt)>0);axis square
subplot(3,2,6)
imagesc(min(abs(output1.M(K,K)),tt)>0);axis square
subplot(3,2,5)
imagesc(min(abs(output1.S(J,J)),tt)>0);axis square
subplot(3,2,6)
imagesc(min(abs(output1.M(J,J)),tt)>0);axis square

figure(41); imagesc(set(J,:)); colormap gray

figure(42)
imagesc(indvalues(J)'); colormap jet;
pbaspect([50 length(industry) 1])

figure(50);
for i=1:size(output1.atoms_u,2)
    clf
    hist(indvalues(output1.atoms_u(:,i)~=0));
    pause
end
% 
% figure(4);clf;
% subplot(2,3,1);
% imagesc(abs(output.S)>1e-2); colormap gray
% axis square;
% subplot(2,3,2);
% imagesc(abs(output.M)>1e-2);
% axis square;
% subplot(2,3,3);
% imagesc(Ac);
% axis square;
% subplot(2,3,4);
% imagesc(abs(output.S(J,J))>1e-2);
% axis square;
% subplot(2,3,5);
% imagesc(abs(output.M(J,J))>1e-2);
% axis square;
% subplot(2,3,6);
% imagesc(Ac(J,J));
% axis square;
% 
% 
% %%
% nb = size(output.atoms_u,2);
% bset = false(p,1);
% for i = 1:nb
%     bset = bset | abs(output.atoms_u(:,i))>1e-10;
% end
% A = Sigma(bset,bset);
% B = Sigma(bset,~bset);
% C = Sigma(~bset,~bset);
% SigBlock = A - B*(inv(C))*B';
% [V,D] = eig(SigBlock);
% ee = diag(D);
% 
% figure(41);clf;
% plot(sort(ee),'.');
% 
% 
% 
