%% EXPERIMENT ON NET1 DREAM5
clear all; clc;
%%
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
addpath('../prox');
addpath('../../glasso-matlab');
addpath ../../ml-100k/

R=importdata('u.data');
IT=importdata('u.item');
GE=importdata('u.genre');

keyboard;
S=cov(F.data);

figure(1);clf;
imagesc(max(abs(S(:)))-abs(S));
colormap bone;

%%
%% our norm psd with decomposition S-M sparse_omega_lgm
po=size(S,1);
p=po;
param.max_nb_main_loop=2;%2;%1000
%%
% param.f=4;
% param.verbose=1;
% inputData.X1=S^.5;
% inputData.X2=S^.5;
%%  
param.f=5;
param.verbose=1;
inputData.X=S;
%%
inputData.Y=-eye(po);

%% reg param
beta=.5;
param.cardfun=((1:p).^beta)/p^beta;
param.cardfun(1)=inf;
param.cardfun(20:end);
lam=.5;
gam=.5;
param.lambda=lam;
param.mu=gam;

%% init with graphical lasso
% rho=1;
% tol=1e-6;
% maxIt=10;
% keyboard;
% [Theta W] = graphicalLasso(S, rho, maxIt, tol);
% keyboard;

%% blocks
[Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
obj_l1_om=hist.obj(end);
save('ml','k','p','n','inputData','Dfull','Dmargo', ...
'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output');

%% reconstruction l1+om
if ~isempty(ActiveSet.alpha)
    Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
    nl=size(ActiveSet.atoms,2);
    Dfin=zeros(p+nl);
    Dfin(1:nl,1:nl)=eye(nl);
    Dfin((nl+1):(nl+p),(nl+1):(nl+p))=-Z1;
    Dfin(1:nl,(nl+1):(nl+p))=Uso';
    Dfin((nl+1):(nl+p),1:nl)=Uso;
else
    Dfin=Z1;
end

keyboard;

%% tr+l1
param.lambda=.6; %lamda ~ 2/k*mu
param.mu=0.1;
param.max_nb_main_loop=2;
param.niterPS=10000;
param.cardfun=inf*ones(1,p);
param.cardfun(p)=1;
[Z_tr Z1_tr Z2_tr ActiveSet_tr hist_tr param_tr flag_tr output_tr] = cgan_l1_omega(inputData,param);

%% reconstruction tr+l1
if ~isempty(ActiveSet_tr.alpha)
    Uso_tr=bsxfun(@times,sqrt(ActiveSet_tr.alpha)',ActiveSet_tr.atoms);
    nl_tr=size(ActiveSet_tr.atoms,2);
    Dfin_tr=zeros(p+nl_tr);
    Dfin_tr(1:nl_tr,1:nl_tr)=eye(nl_tr);
    Dfin_tr((nl_tr+1):(nl_tr+p),(nl_tr+1):(nl_tr+p))=-Z1_tr;
    Dfin_tr(1:nl_tr,(nl_tr+1):(nl_tr+p))=Uso_tr';
    Dfin_tr((nl_tr+1):(nl_tr+p),1:nl_tr)=Uso_tr;
else
    Dfin_tr=Z1_tr;
end



save('ml','k','p','n','inputData','Dfull','Dmargo', ...
'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output',...
'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr');

obj_l1_om
hist_tr.obj(end)