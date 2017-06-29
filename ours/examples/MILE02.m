%% EXPERIMENT ON NET1 DREAM5
clear all; clc;
%%
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
addpath('../prox');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');
addpath('../../glasso-matlab');

run('pp_MILE.m');

p=size(S,1);

%%
param.f=4;
param.verbose=1;
inputData.X1=S^.5;
inputData.X2=S^.5;
%%  
% param.f=5;
% param.verbose=1;
% inputData.X=S;
% %%

%% our norm psd with decomposition S-M sparse_omega_lgm
inputData.Y=-eye(p);
param.cardfun=inf*ones(1,p);
param.cardfun(100)=1;
lam=.11;
gam=.005;
param.lambda=lam;
param.mu=gam;
param.max_nb_main_loop=1;


%% blocks
% [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
[Z L S D ActiveSet hist param flag output] = cgan_l1_omega_02(inputData,param);
Z1=S+D;
Z2=L;

%% objective
figure(1);clf;
plot(hist.obj0);

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


%%

figure(10);clf;
subplot(1,4,1);
imagesc(abs(Dfin)>1e-10);
axis square;
subplot(1,4,2);
imagesc(abs(Dfin));
axis square;
subplot(1,4,3);
imagesc(abs(Z1));
title('S')
axis square;
subplot(1,4,4);
imagesc(abs(Z2));
title('L')
axis square;
colormap hot;

%%
%% reorder

[I]=grayorder(Uso~=0);
UsoI=Uso(I,:);
Z2II=Z2(I,I);
Z1II=Z1(I,I);


figure(11);clf;
subplot(2,3,1);
imagesc(abs(UsoI));
axis square;
subplot(2,3,2);
imagesc(abs(Z1II));
title('L')
axis square;
subplot(2,3,3);
imagesc(abs(Z2II));
axis square;
subplot(2,3,4);
imagesc(abs(UsoI)>1e-10);
axis square;
subplot(2,3,5);
imagesc(abs(Z1II)>1e-10);
title('L')
axis square;
subplot(2,3,6);
imagesc(abs(Z2II)>1e-10);
colormap hot
axis square;


% save('mile','Dfin', ...
% 'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output');

% %% tr+l1
% param.lambda=.6; %lamda ~ 2/k*mu
% param.mu=0.1;
% param.max_nb_main_loop=2;
% param.niterPS=10000;
% param.cardfun=inf*ones(1,p);
% param.cardfun(p)=1;
% [Z_tr Z1_tr Z2_tr ActiveSet_tr hist_tr param_tr flag_tr output_tr] = cgan_l1_omega(inputData,param);
% 
% 
% 
% 
% %% reconstruction tr+l1
% if ~isempty(ActiveSet_tr.alpha)
%     Uso_tr=bsxfun(@times,sqrt(ActiveSet_tr.alpha)',ActiveSet_tr.atoms);
%     nl_tr=size(ActiveSet_tr.atoms,2);
%     Dfin_tr=zeros(p+nl_tr);
%     Dfin_tr(1:nl_tr,1:nl_tr)=eye(nl_tr);
%     Dfin_tr((nl_tr+1):(nl_tr+p),(nl_tr+1):(nl_tr+p))=-Z1_tr;
%     Dfin_tr(1:nl_tr,(nl_tr+1):(nl_tr+p))=Uso_tr';
%     Dfin_tr((nl_tr+1):(nl_tr+p),1:nl_tr)=Uso_tr;
% else
%     Dfin_tr=Z1_tr;
% end
% 
% obj_l1_om
% hist_tr.obj(end)