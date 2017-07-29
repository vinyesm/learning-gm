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

addpath('../../spca_am-master');

HOME = '../../code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

fig=0;

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
fig=0;
%% OURS
param.mu=.1;
param.lambda=.1;

%% LVGGM Chandrasekaran S-L, (Sparse-Low Rank)

%%
%% set up SDP data in SDPT3 format
%%
tic
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
time_chandra=toc;

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

if fig
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
end

ActiveSet.I={};
ActiveSet.k={};
param.k=100;
param.verbose=1;
param.max_nb_main_loop=100;
param.cardfun=inf*ones(1,p);
param.cardfun(param.k)=1;
tic
[M,S,L,U,hist,ActiveSet] = logdetPPA_l1_omega(Sigma,param,ActiveSet);
time_solver=toc

if fig
    figure(2);clf;
    plot(hist.obj_sup);
end

%
if ~isempty(ActiveSet.I)
    Uso=[];
    for i=1:length(ActiveSet.I);
        Uso=[Uso U{i}];
    end
    nl=size(Uso,2);
    Dfin2=zeros(p+nl);
    Dfin2(1:nl,1:nl)=eye(nl);
    Dfin2((nl+1):(nl+p),(nl+1):(nl+p))=S;
    Dfin2(1:nl,(nl+1):(nl+p))=Uso';
    Dfin2((nl+1):(nl+p),1:nl)=Uso;
else
    Dfin2=S;
end

if fig
    figure(3);clf;
    subplot(1,2,1)
    imagesc(abs(Dfin2));
    pbaspect([1 1 1]);
    title('estimated complete conc. mat.');
    colorbar
    subplot(1,2,2)
    imagesc(abs(Dfin2)>1e-15);
    pbaspect([1 1 1]);
    title('estimated support');
    colorbar
end

