%
% Chandrasekaran on Breast Cancer Data


%% add paths
clc; clear all; close all;

 HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
% HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

%% data
load('BC.mat')

%% LVGGM Chandrasekaran S-L, (Sparse-Low Rank)

%% parama
lambda=.05;
mu=.001;

%% set up SDP data in SDPT3 format
%%
tic
n=size(Sigma,1);
rho=mu;
invD = speye(n,n); 

%
n2 = n*(n+1)/2;
b = zeros(n2,1);
C{1} = Sigma;
blk{1,1} = 's'; blk{1,2} = n;
[Iall,Jall] = find(triu(ones(n)));
tmp = [Iall,Jall];
m2 = size(tmp,1);
Icomp = tmp(:,1); Jcomp = tmp(:,2);
Itmp  = Icomp + Jcomp.*(Jcomp-1)/2;
Atmp  = spconvert([Itmp,[1:m2]',ones(m2,1); n2,m2,0]);
At{1} = Atmp;
%%
beta = lambda;
blk{2,1} = 's'; blk{2,2} = n;
At{2,1}  = Atmp;
C{2,1}   = beta*speye(n,n);
%%
blk{3,1} = 'l'; blk{3,2} = 2*n2;
Identity = speye(n2);
At{3,1} = [-Identity,Identity]';
idx = find(Icomp == Jcomp);
ee  = sqrt(2)*ones(m2,1);
if ~isempty(idx); ee(idx) = ones(length(idx),1); end
C{3,1} = rho*[ee; ee];
%fprintf('\n Set up data time = %3.2f',etime(clock,ttime));
runPPA = 1;
if (runPPA)
    OPTIONS.smoothing  = 1;
    OPTIONS.scale_data = 0; %% or 2;
    OPTIONS.plotyes    = 0;
    OPTIONS.tol        = 1e-12;
    mu = [1; 0; 0];
    [obj,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,mu,OPTIONS);
    obj = sum(sum(Sigma.*X{1}))-sum(log(eig(X{1})))+rho*sum(sum(abs(X{1})));
    X1 = invD*X{1}*invD; X1 = 0.5*(X1+X1');
    X2 = invD*X{2}*invD; X2 = 0.5*(X2+X2');
end
time_chandra=toc;

%% solution
Ssl=X1+X2;
Lsl=X2;

%%
figure(1);
subplot(1,2,1)
imagesc(abs(Ssl-diag(diag(Ssl)))); colormap hot;
axis square
subplot(1,2,2)
imagesc(abs(Lsl)); colormap hot;
axis square
