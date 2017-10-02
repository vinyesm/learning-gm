%figures of the experiment output
%
addpath ../reorder/
load('../genedata/BC.mat')
load('exp_50_100_100.mat')


%K1=S1-M1;
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
