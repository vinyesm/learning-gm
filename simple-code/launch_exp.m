% clear all;close all;clc
% 
% addpath ../ours/TPower_1.0/misc/
% addpath ../ours/TPower_1.0/algorithms/TPower/
% addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
% addpath ../spams-matlab-v2.6/build/
% 
% %HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
%  HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
% addpath(strcat(HOME,'/solver/'))
% addpath(strcat(HOME,'/solver/mexfun'))
% addpath(strcat(HOME,'/util/'))
% load('../genedata/BC.mat')
% ks =  [25 50 100]
% las = [.1 .1  .01  .01  .01];
% mus = [.1 .01 .01  .005 .001];

function [] = launch_exp(lambda, mu, k, Sigma_train, Sigma_test, savename)
    param.k=k;
    param.epsObj=1e-5;
    param.lambda=lambda;
    param.mu=mu;
    param.maxIter=20;
    param.maxNbBlocks=100;
    param.verbose=2;
    %% logdetOmegaL1 initialised with true support
    [S1,M1,L1,U1,hist_ch1,set1] = logdetOmegaL1(Sigma_train,param,inf,Sigma_test);
    %loglikelihood = @(K,Sig) log(det(K))-trace(K*Sig);
    K = S1-M1;
    %nb param of t
    logltrain = log(det(K))-trace(K*Sigma_train);
    logltest = log(det(K))-trace(K*Sigma_test);
    fprintf('loglikelihood train = %f       loglikelihood test = % f\n', logltrain, logltest )
    save(['exp_' savename], 'S1','M1','L1','U1','hist_ch1','set1', 'param','logltrain','logltest')
end




