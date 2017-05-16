%% EXPERIMENT ON NET1 DREAM5
clear all; clc;
%%
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
addpath('../prox');
addpath('../../glasso-matlab');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

run('pp_SP100.m');
% run('pp_movielens_2.m');
% k=20;
%%
%% our norm psd with decomposition S-M sparse_omega_lgm
po=size(S,1);
p=po;
param.max_nb_main_loop=50;%2;%1000
%%
param.f=4;
param.verbose=1;
inputData.X1=real(S^.5);
inputData.X2=real(S^.5);

%%  
% param.f=5;
% param.verbose=1;
% inputData.X=S;
%%
inputData.Y=-eye(po);

%% reg param

% param.cardfun=inf*ones(1,p);
% param.cardfun(k)=1;
beta=.5;
param.cardfun=((1:p).^beta)/p^beta;
param.cardfun(1)=inf;
param.cardfun(50:end)=inf;
param.lambda=.2;
param.mu=.2;


%% blocks
% [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
% 
% save(['SP100_' num2str(param.f) '_' num2str(param.mu) '_' num2str(param.lambda) ],'param', ...
% 'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output');


%% tr+l1
%% Starting solution
param.cardfun=inf*ones(1,p);
param.cardfun(p)=1;

param.lambda=.008; %lamda ~ 2/k*mu
param.mu=.003;
param.max_nb_main_loop=2;
param.niterPS=10000;
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

figure(3);clf;
subplot(1,4,1)
imagesc(abs(Dfin_tr));
pbaspect([1 1 1]);
subplot(1,4,2)
imagesc(abs(Dfin_tr)>0);
pbaspect([1 1 1]);
subplot(1,4,3)
imagesc(abs(Z_tr));
pbaspect([1 1 1]);
subplot(1,4,4)
imagesc(abs(Z2_tr));
pbaspect([1 1 1]);

save(['SP100_' num2str(param.f) '_' num2str(param.mu) '_' num2str(param.lambda) ],'param', ...
'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr');
% save('ml','k','p','n','inputData','Dfull','Dmargo', ...
% 'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output',...
% 'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr');

% obj_l1_om
% hist_tr.obj(end)