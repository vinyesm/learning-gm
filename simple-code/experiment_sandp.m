clear all;close all;clc

% diviser les donnees

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/

%HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
% HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
HOME = '../code-from-Kim-Chuan/LogdetPPA-0';
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

load('../sandp500/sandp500.mat')
addpath ../reorder/



%ks =  [size(Sigma_train,1) 150];
ks = [size(Sigma_train,1) 100];
las = .7;
mus = .2;
% las = [5   2   2  1  1 .1];
% mus = [.5 .5 .1 .5  .1 .05];


% parfor i=1:length(las)
for i=1:length(las)
    %chandra
    for k=1:length(ks)
        k0 = ks(k);
        la0 = las(i);
        mu0 = mus(i);
        savename = [num2str(k0) '_' num2str(la0*1000) '_' num2str(mu0*1000)];
        % nombre de parametres, flag si finished et logvraissemblace sur
        % les autres donnees
        
        launch_exp(la0, mu0, k0, Sigma_train, Sigma_test, [savename '_train'])
        launch_exp(la0, mu0, k0, Sigma_test, Sigma_train, [savename '_test'])
    end
end

%%
% J=grayorder(full(set1~=0));
% figure(30);clf;
% U=[U1{:}];
% imagesc(U(J,:))
% 
% figure(31);
% imagesc(indvalues(J)');
% 
% industry2 =industry(I);
% indvalues2=indvalues(I);
% [val,ind]=sort(industry2);
% 
% figure(32);clf;
% U=[U1{:}];
% imagesc([U(ind,:) indvalues2(ind)'/20])
% 
% % confusion matrix
% Iind = double(bsxfun(@eq, repmat(indvalues2',1,11), 1:11));
% Jaccard = (set1'*Iind)./(123-(1-set1)'*(1-Iind));
% conf = (set1'*Iind)./repmat(sum(Iind),5,1); % recall
% conf2 = (set1'*Iind)./repmat(sum(set1)',1,11); 
% 
% figure(33);
% imagesc(conf); colormap gray
% 
% figure(34);
% imagesc(conf2); colormap gray

%%
A=load('exp_100_700_200_train.mat');
B=load('exp_100_700_200_test.mat');
C=load('exp_123_700_200_train.mat');
D=load('exp_123_700_200_test.mat');

full(A.set1'*B.set1)
esperance = (100/123)^2*123;

UA = normc([A.U1{:}]);
[UAA, ~, ~]= svd(UA);
UB = normc([B.U1{:}]);
[UBB, ~, ~]= svd(UB);
[U0, S0, V0]= svds(UAA'*UBB);

J=grayorder(full(A.set1~=0));
figure(40);
subplot(1,2,1);
imagesc(A.set1(J,:)); colormap gray
subplot(1,2,2);
imagesc(B.set1(J,:)); colormap gray

norm(C.M1-D.M1,'fro')^2/(norm(C.M1,'fro')*norm(D.M1,'fro'));
UC=normc(C.U1{1});
UD=normc(D.U1{1});
[U, S, V]= svds(UC'*UD);
