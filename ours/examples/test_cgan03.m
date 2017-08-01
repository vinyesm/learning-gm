
% TESTING NEW FUNCTION cgan_l1_omega.m
% EXPERIMENT ON  LARGE BLOCKS AND SPARSE MATRIX ON OBSRVED VAR
%non oriented graph structure learning with latent variables
% 
% We build the complete model and sample from it
% We assume latept variables ndependept

%% add paths
% clc
% clc; clear all; close all;
% addpath('../main');
% addpath('../active-set');
% addpath('../atom-selection');
% addpath('../utils');
% addpath('../other');
% addpath('../prox');
% addpath('../TPower_1.0');
% addpath('../TPower_1.0/algorithms/TPower/');
% addpath('../TPower_1.0/misc/');
% addpath('../../spams-matlab-v2.6/build/');
% 
% % %% data
% run('../../toy-data/toy_overlap.m');k=15;

%% our norm psd with decomposition S-M sparse_omega_lgm
p=po;
param.verbose=1;
inputData.Y=-eye(po);

%%
param.f=4;
inputData.X1=S^.5;
inputData.X2=S^.5;

%%
inputData.Y=-eye(po);
param.cardfun=inf*ones(1,p);
param.cardfun(k)=1;

%%% n=5000
param.lambda=.7; %lamda ~ 2/k*mu
param.mu=.1;

param.sloppy=0;
param.max_nb_main_loop=1500;
param.niterPS=100;


%%
% %% Starting solution

ActiveSet.max_atom_count_reached=0;
ActiveSet.I={};
ActiveSet.alpha= [];
ActiveSet.atoms=pl;
ActiveSet.atom_count = pl;
ActiveSet.k=mat2cell(ks,1,ones(1,length(ks)));
ActiveSet.alpha=sum(Dol.^2)';
ActiveSet.atoms=sparse(bsxfun(@rdivide, Dol, sqrt(sum(Dol.^2))));
for i=1:pl
    ActiveSet.I{i}=find(Dol(:,i));
end

cf=inf*ones(1,length(ks));
for j=1:length(ks)
    cf(j)=min(param.cardfun(ks(j):end));
end

ActiveSet.atoms=bsxfun(@rdivide,ActiveSet.atoms(:,1:ActiveSet.atom_count),sqrt(cf));
ActiveSet.alpha=ActiveSet.alpha.*cf';

I=speye(p)==1;
startingZ.Z1=-Doo;
startingZ.Z2=Dol*Dol';



%% blocks
%[Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param,startingZ,ActiveSet);
% obj_l1_om=hist.obj(end);
% keyboard;

%% tr+l1
% [Zb Lb Sb ActiveSet2 hist2 param2 flag2 output2] = cgan_l1_omega_03(inputData,param,startingZ,ActiveSet);
[Zb Lb Sb ActiveSet2 hist2 param2 flag2 output2] = cgan_l1_omega_03(inputData,param);%,startingZ,ActiveSet);
Z1b=Sb;
Z2b=Lb;

%% objective
figure(1);clf;
semilogy(hist2.obj0);
title('objective')

figure(2);clf;
semilogy(hist2.dg);
title('dg sub-problems PS')

figure(3);clf;
semilogy(hist2.dg_S,'r');hold on;
semilogy(hist2.dg_L,'b');hold on;
legend({'S dg' 'L dg'});
title('duality gaps of subproblems')

figure(4);clf;
semilogy(hist2.dg_global);
title('dg global')


%% LOAD HERE lggm.mat

%% reconstruction l1+om
% if ~isempty(ActiveSet.alpha)
%     Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
%     nl=size(ActiveSet.atoms,2);
%     Dfin=zeros(p+nl);
%     Dfin(1:nl,1:nl)=eye(nl);
%     Dfin((nl+1):(nl+p),(nl+1):(nl+p))=-Z1;
%     Dfin(1:nl,(nl+1):(nl+p))=Uso';
%     Dfin((nl+1):(nl+p),1:nl)=Uso;
% else
%     Dfin=Z1;
% end

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

%% reconstruction l1+om
if ~isempty(ActiveSet2.alpha)
    Uso2=bsxfun(@times,sqrt(ActiveSet2.alpha)',ActiveSet2.atoms);
    nl=size(ActiveSet2.atoms,2);
    Dfin2=zeros(p+nl);
    Dfin2(1:nl,1:nl)=eye(nl);
    Dfin2((nl+1):(nl+p),(nl+1):(nl+p))=-Z1b;
    Dfin2(1:nl,(nl+1):(nl+p))=Uso2';
    Dfin2((nl+1):(nl+p),1:nl)=Uso2;
else
    Dfin2=Z1b;
end


%%

figure(5);clf;
subplot(2,2,1)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
colorbar
subplot(2,2,2)
imagesc(abs(Dfin2));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(2,2,3)
imagesc(abs(Dfull)>1e-15);
pbaspect([1 1 1]);
title('true support');
colorbar
subplot(2,2,4)
imagesc(abs(Dfin2)>1e-15);
pbaspect([1 1 1]);
title('estimated support');
colorbar


% %%
% figure(3);clf;
% subplot(3,2,1)
% axis off;
% text(0,.5,descr_gen)
% pbaspect([1 1 1]);
% subplot(3,2,2)
% axis off;
% text(0,.5,descr_par)
% pbaspect([1 1 1]);
% subplot(3,2,3)
% imagesc(abs(Dfull));
% pbaspect([1 1 1]);
% title('true complete conc. mat.');
% colorbar
% subplot(3,2,4)
% imagesc(abs(Dfin));
% pbaspect([1 1 1]);
% title('estimated complete conc. mat.');
% colorbar
% subplot(3,2,5)
% imagesc(abs(Dfull)>1e-15);
% pbaspect([1 1 1]);
% title('true support');
% colorbar
% subplot(3,2,6)
% imagesc(abs(Dfin)>1e-15);
% pbaspect([1 1 1]);
% title('estimated support');
% colorbar
% 
% THRESH=1e-3;
% 
% figure(6);clf;
% subplot(3,2,1)
% imagesc(abs(Dfull));
% pbaspect([1 1 1]);
% title('true complete conc. mat.');
% colorbar
% subplot(3,2,2)
% imagesc(abs(Dfull));
% pbaspect([1 1 1]);
% title('true support');
% colorbar
% subplot(3,2,3)
% imagesc(abs(Dfin).*(abs(Dfin)>THRESH));
% pbaspect([1 1 1]);
% title('estimated complete conc. mat.');
% colorbar
% subplot(3,2,4)
% imagesc(abs(Dfin)>THRESH);
% pbaspect([1 1 1]);
% title('estimated support');
% colorbar
% subplot(3,2,5)
% imagesc(abs(Dfin_tr).*(abs(Dfin_tr)>THRESH));
% pbaspect([1 1 1]);
% title('estimated complete conc. mat.');
% colorbar
% subplot(3,2,6)
% imagesc(abs(Dfin_tr)>THRESH);
% pbaspect([1 1 1]);
% title('estimated support');
% colorbar
% colormap jet
% 
% % myColorMap = hot; % Make a copy of jet.
% % % Assign white (all 1's) to black (the first row in myColorMap).
% % myColorMap(1, :) = [1 1 1];
% % colormap(myColorMap); % Apply the colormap 
% 
% 
% % nb good edges
% Z1T=Z1.*(abs(Z1)>THRESH);
% Z1trT=Z1_tr.*(abs(Z1_tr)>THRESH);
% 
% rec1=1-sum(sum(sign(Doo)+sign(Z1T)))/p^2;
% rec1_tr=1-sum(sum(sign(Doo)+sign(Z1trT)))/p^2;
% fprintf('rec1=%f  rec1_tr=%f\n',rec1,rec1_tr);
% 
% save('TOY_OVERLAP2', 'Dmargo', 'Doo', 'Dfull', 'S', ...
%                      'Z', 'Z1', 'Z2', 'ActiveSet', 'hist', 'param', 'flag', 'output', 'Dfin', ...
%                      'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr', 'Dfin_tr', ...
%                      'THRESH');
%                      