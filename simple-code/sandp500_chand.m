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
load('../sandp500/sandp500.mat')

%Sigma = genDataGraph(Ac,50000);

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


%% logdetOmegaL1 initialised with true support
% lambda = 10;
% mus = .02;
lambda = 0.01;
mus = 0.0005;
%mus = 2.^linspace(-16,1,30);
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
end


%%
figure(1);clf;
subplot(1,2,1);
imagesc(abs(Ssl)>1e-2);
axis square;
subplot(1,2,2);
imagesc(abs(Lsl)>1e-2);colormap gray; 
axis square;

%%
% Z = linkage(Lsl,'ward');
% % [H,T,OUTPERM] = dendrogram(Z) ;
% %[Cres,I]=order_of_tree(Z);
% I = optimalleaforder(Z,pdist(Lsl));

% Acc = Ac + eye(size(Ac,1));
% Z = linkage(full(Acc),'ward');
% % [H,T,OUTPERM] = dendrogram(Z) ;
% %[Cres,I]=order_of_tree(Z);
% I = optimalleaforder(Z,pdist(full(Acc)));
% figure(31);
% imagesc(Ac(I,I));colormap gray; 



figure(2);clf;
subplot(1,2,1);
imagesc(min(abs(Ssl(I,I)),10));colormap jet; 
axis square;
subplot(1,2,2);
imagesc(min(abs(Lsl(I,I)),10));
axis square;


%%


figure(3);clf;
subplot(1,2,1);
imagesc(min(abs(Ssl),10));
axis square;
subplot(1,2,2);
imagesc(min(abs(Lsl),10));colormap jet; 
axis square;


[V,D]=eig(Lsl);
