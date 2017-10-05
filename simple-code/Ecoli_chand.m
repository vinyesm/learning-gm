clear all;close all;clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/

HOME = '../code-from-Kim-Chuan/LogdetPPA-0';
%HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
% HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

addpath ../reorder/

%load('../genedata/BC.mat')
%run('../DREAM5/Dream/pp_net1');
load('../DREAM5/Dream/pp_net1.mat')

%%
n=size(Sigma,1);
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


OPTIONS.smoothing  = 1;
OPTIONS.scale_data = 0; %% or 2;
OPTIONS.plotyes    = 0;
OPTIONS.tol        = 1e-12;
eta = [1; 0; 0];



%%

pos = sum(sum((Ac==1)));
neg = sum(sum((Ac==0)));
tpr = [];
fpr = [];

%% logdetOmegaL1 initialised with true support
lambda = .1;
mus = 2.^linspace(-16,1,30);
for mu=mus
    %%
    blk{2,1} = 's'; blk{2,2} = n;
    At{2,1}  = Atmp;
    C{2,1}   = lambda*speye(n,n);
    %%
    blk{3,1} = 'l'; blk{3,2} = 2*n2;
    Identity = speye(n2);
    At{3,1} = [-Identity,Identity]';
    idx = find(Icomp == Jcomp);
    ee  = sqrt(2)*ones(m2,1);
    if ~isempty(idx); ee(idx) = ones(length(idx),1); end
    C{3,1} = mu*[ee; ee];
    [obj,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,eta,OPTIONS);
    %obj = sum(sum(Sigma.*X{1}))-sum(log(eig(X{1})))+mu*sum(sum(abs(X{1})));
    X1 = invD*X{1}*invD; X1 = 0.5*(X1+X1');
    X2 = invD*X{2}*invD; X2 = 0.5*(X2+X2');
    %% solution
    Ssl=X1+X2;
    Lsl=X2;

    Sc = abs(Ssl-diag(diag(Ssl)))>0;
    tp = sum(sum((Ac==1) & Sc));
    fp = sum(sum((Ac==0) & Sc));
    tpr = [tpr, tp/pos];
    fpr = [fpr, fp/(fp+neg)];
end


figure(1);clf
plot(fpr,tpr,'.');hold on;
plot([0 1],[0,1],'r-'); 
axis([0 1 0 1])
