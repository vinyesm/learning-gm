load('../DREAM5/Dream/pp_net1.mat')
% addpath ../ours/examples/
% addpath ../ml-100k/
% addpath ../ours/utils/

% run('pp_movielens.m');
% Sigma=S;

param.k=20;
param.epsObj=1e-16;
% param.lambda=.03;
% param.mu=.02;
param.lambda=300;
param.mu=20;
param.maxIter=500;
param.maxNbAtoms=1000;
param.verbose=1;
param.epsStop=1e-4;

p = size(Sigma,1);
inputData.X = Sigma^0.5;
inputData.Y = eye(p);


[ output, hist ] = regOmegaL1( inputData, param, inf );

%%
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
imagesc(output.M);
% figure(18)
% imagesc(M);

figure(19)
imagesc(-full(output.S));
% figure(20)
% imagesc(S);

%%
ua = unique(abs(output.atoms_u')>1e-10, 'rows');
set = ua';
figure(40); imagesc(set); colormap gray

% %%
% [J]=grayorder(full(set~=0));
% figure(3);clf;
% subplot(2,2,1);
% imagesc(min(abs(output.S),10));
% axis square;
% subplot(2,2,2);
% imagesc(min(abs(output.M),10));
% axis square;
% subplot(2,2,3);
% imagesc(min(abs(output.S(J,J)),10));
% axis square;
% subplot(2,2,4);
% imagesc(min(abs(output.M(J,J)),10));
% axis square;
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