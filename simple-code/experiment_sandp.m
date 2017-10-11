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



%ks =  [size(Sigma_train,1) 150];
ks = [size(Sigma_train,1) 150 50 100];
las = 2;
mus = .5;
%las = [5   2   2  1  1 .1];
%mus = [.5 .5 .1 .5  .1 .05];


parfor i=1:length(las)
    %chandra
    for k=1:length(ks)
        k0 = ks(k);
        la0 = las(i);
        mu0 = mus(i)
        savename = [num2str(k0) '_' num2str(la0*1000) '_' num2str(mu0*1000)];
        % nombre de parametres, flag si finished et logvraissemblace sur
        % les autres donnees
        
        launch_exp(la0, mu0, k0, Sigma_train, Sigma_test, savename)
    end
end