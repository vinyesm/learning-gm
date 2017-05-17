%
% EXPERIMEpt ON non oriepted graph structure learning with latept variables
%
% We build the complete model and sample from it
% We assume latept variables ndependept

%% add paths
clc; clear all; close all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
% addpath('../other');
addpath('../prox');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

% HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

% addpath LogdetPPA-0/;
% addpath LogdetPPA-0/solver/;
% addpath LogdetPPA-0/solver/mexfun/;
% addpath LogdetPPA-0/util/;


%% data
% run('pp_movielens.m');
run('pp_MILE.m');
% run('pp_MILE_bis.m');
% run('pp_movielens_2.m');

% WHEN 500 genes
% lambda=2.5;
% mu=.02;

% WHEN 200 genes
lambda=1;
mu=.02;

%% LVGGM Chandrasekaran S-L, (Sparse-Low Rank)

%%
%% set up SDP data in SDPT3 format
%%
Sigma=S;
n=size(S,1);
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

%% solution
Ssl=X1+X2;
Lsl=X2;

%Usl=chol(Lsl);
[VV DD]=eig(Lsl);
dd=diag(DD);
VV=VV(:,dd>1e-8);
dd=dd(dd>1e-8);
Usl=VV*diag(sqrt(dd));

p=n;
nl=length(dd);

Dsl=zeros(p+nl);
Dsl(1:nl,1:nl)=eye(nl);
Dsl((nl+1):(nl+p),(nl+1):(nl+p))=Ssl;
Dsl(1:nl,(nl+1):(nl+p))=Usl';
Dsl((nl+1):(nl+p),1:nl)=Usl;


figure(3);clf;
subplot(1,4,1);
imagesc(min(abs(Dsl),1));
title('Complete');
colormap hot
axis square
subplot(1,4,2);
imagesc(abs(Dsl)>1e-3);
title('support Complete');
axis square
subplot(1,4,3);
imagesc(abs(Ssl));
title('S');
axis square
subplot(1,4,4);
imagesc(abs(Lsl));colormap hot
title('L');
axis square

save('Lsl','Lsl', 'Ssl');


%%%___________________________
%%
% S2=inv(2*eigs(Lsl,1, 'la')*eye(p)-Lsl);
p=size(Lsl,1);
param.verbose=1;
%
param.f=1; %prox
param.verbose=1;
inputData.Y=eye(p)+Lsl;

%
% param.f=4;
% param.verbose=1;
% inputData.X1=S2^.5;
% inputData.X2=S2^.5;
% inputData.Y=-eye(p);

%
% param.f=5;
% inputData.X=S;
% inputData.Y=-eye(p);

%
% reg param
% beta=.5;
% param.cardfun=((1:p).^beta)/p^beta;
% param.cardfun(1:5)=inf;
% param.cardfun(150:end)=inf;
param.cardfun=inf*ones(1,p);
param.cardfun(100)=1;
lam=10;
% lam=12;
gam=100;
param.lambda=lam;
param.mu=gam;
param.max_nb_main_loop=100;

%% blocks
[Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);


%% 
if ~isempty(ActiveSet.alpha)
    Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
else
    Uso=zeros(p,1);
    Dfin=Z1;
end

%% reorder

[I]=grayorder(Uso~=0);

% Z = linkage(full(Uso),'ward');
% [Cres,I]=order_of_tree(Z);
% 
genesI0=genesI;
genesI=genesI(I);
UsoI=Uso(I,:);
Z2II=Z2(I,I);
imp_idxII=imp_idxI(I);
%%

figure(1);clf;
subplot(2,3,1);
imagesc(abs(Uso)>1e-10);
axis square;
subplot(2,3,2);
imagesc(abs(Z2));
title('L')
axis square;
subplot(2,3,3);
imagesc(imp_idxI');
axis square;
subplot(2,3,4);
imagesc(abs(UsoI)>1e-10);
axis square;
subplot(2,3,5);
imagesc(abs(Z2II));
title('L')
axis square;
subplot(2,3,6);
imagesc(imp_idxII');
colormap hot
axis square;

figure(4);clf;
subplot(2,3,1);
imagesc(abs(Ssl)>1e-3);
title('S');
axis square
subplot(2,3,2);
imagesc(abs(Lsl));colormap hot
title('L');
axis square
subplot(2,3,3);
imagesc(imp_idxI');
axis square;
subplot(2,3,4);
imagesc(abs(Ssl(I,I))>1e-3);
title('S');
axis square
subplot(2,3,5);
imagesc(abs(Lsl(I,I)));colormap hot
title('L');
axis square
subplot(2,3,6);
imagesc(imp_idxII');
colormap hot
axis square;





