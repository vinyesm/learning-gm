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
load('../genedata/BC.mat')

ks =  [25 50 100];
las = [.1 .1  .01  .01  .01];
mus = [.1 .01 .01  .005 .001];

parfor i=1:length(las)
    for k=1:length(ks)
        k0 = ks(k);
        la0 = las(i);
        mu0 = mus(i)
        savename = [num2str(k0) '_' num2str(la0*1000) '_' num2str(mu0*1000)];
        launch_exp(la0, mu0, k0, Sigma, savename)
    end
end