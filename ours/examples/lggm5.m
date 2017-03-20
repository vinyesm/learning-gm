% TESTING NEW FUNCTION cgan_l1_omega.m
% EXPERIMENT ON BLOCKS OF DIFFERENT SIZES
%non oriented graph structure learning with latent variables
% 
% We build the complete model and sample from it
% We assume latept variables ndependept

%% add paths
clc; clear all; close all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
addpath('../other');
addpath('../prox');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

%% data
run('../../toy-data/three_sizes_blocks_20x10x5.m');

%% our norm psd with decomposition S-M sparse_omega_lgm
p=po;
param.f=4;
param.diag=0;
param.PSD=true;
param.max_nb_main_loop=100;%2;%1000
param.powerIter=500;
param.stPtPowerIter=1000;
param.niterPS=5000;%10000;%5000
param.epsStop=1e-8;
param.PSdualityEpsilon=1e-8;
param.k=0;
param.PSmu=0; %strong convexity
param.verbose=1;
param.debug=0;
param.sloppy=0;
param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
inputData.X1=S^.5;
inputData.X2=S^.5;
inputData.Y=-eye(po);
% param.cardfun=inf*ones(1,p);
% param.cardfun(k)=1;
% param.cardfun(p)=1;
beta=0.2;
param.cardfun=((1:p).^beta)./(p^beta);
param.cardfun(1)=inf;

aa=0.01;
param.lambda=aa; %lamda ~ 2/k*mu
param.mu=0.1;


%% blocks
[Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);

%% tr+l1
param.max_nb_main_loop=2;
param.niterPS=1;%5000;
param.cardfun=inf*ones(1,p);
param.cardfun(p)=1;
[Z_tr Z1_tr Z2_tr ActiveSet_tr hist_tr param_tr flag_tr output_tr] = cgan_l1_omega(inputData,param);

%% LOAD HERE lggm.mat

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

%%
descr_gen = {'Formulation Z1+Z2 where ';
    'with sparse + Omega symmetric regularization';
    'solved with our algorithm';
    'unique atomic norm';
    'no psd constraint Z1+Z2';
    ['number of samples n=' num2str(n)];
    };

descr_par = {['\lambda = ' num2str(param.lambda) '  (omega reg)'];
    ['\mu = ' num2str(param.mu) '  (l_1 reg)'];};

%%
figure(3);clf;
subplot(3,2,1)
axis off;
text(0,.5,descr_gen)
pbaspect([1 1 1]);
subplot(3,2,2)
axis off;
text(0,.5,descr_par)
pbaspect([1 1 1]);
subplot(3,2,3)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
colorbar
subplot(3,2,4)
imagesc(abs(Dfin));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(3,2,5)
imagesc(abs(Dfull)>1e-15);
pbaspect([1 1 1]);
title('true support');
colorbar
subplot(3,2,6)
imagesc(abs(Dfin)>1e-15);
pbaspect([1 1 1]);
title('estimated support');
colorbar
%%

figure(4);clf;
subplot(2,2,1);
imagesc(abs(Dmargo));
pbaspect([1 1 1]);
title('true marginal conc. mat.');
colorbar
subplot(2,2,2);
imagesc(abs(Z1+Z2));
pbaspect([1 1 1]);
title('observed conc. mat.');
colorbar
subplot(2,2,3)
imagesc(abs(Dmargo)>1e-15);
pbaspect([1 1 1]);
title('true support');
colorbar
subplot(2,2,4)
imagesc(abs(Z1+Z2)>1e-15);
pbaspect([1 1 1]);
title('estimated support');
colorbar

figure(5);clf
loglog(hist.time,hist.dg,'-','LineWidth',2,'Color',[1 0 0],'DisplayName','dg');hold on;
loglog(hist.time_sup,hist.dg_sup,'-','LineWidth',2,'Color',[0 0 0],'DisplayName','dg sup');hold on;
legend('show','Location','southwest');
grid on
hold off

figure(6);clf;
subplot(3,2,1)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
colorbar
subplot(3,2,2)
imagesc(abs(Dfull)>1e-15);
pbaspect([1 1 1]);
title('true support');
colorbar
subplot(3,2,3)
imagesc(abs(Dfin));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(3,2,4)
imagesc(abs(Dfin)>1e-15);
pbaspect([1 1 1]);
title('estimated support');
colorbar
subplot(3,2,5)
imagesc(abs(Dfin_tr));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(3,2,6)
imagesc(abs(Dfin_tr)>1e-15);
pbaspect([1 1 1]);
title('estimated support');
colorbar


%% saving
%filename = ['lggm2_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.ps'];
filename = ['lggm4_' datestr(clock) '.ps'];
%print ( '-dpsc2', filename, '-f1' )
print ( '-dpsc2', filename, '-append', '-f1' )
print ( '-dpsc2', filename, '-append', '-f2' )
print ( '-dpsc2', filename, '-append', '-f3' )
print ( '-dpsc2', filename, '-append', '-f4' )
print ( '-dpsc2', filename, '-append', '-f5' )

% keyboard
save('lggm4_03_16','p','n','inputData','Dfull','Dmargo', ...
'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output',...
'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr');