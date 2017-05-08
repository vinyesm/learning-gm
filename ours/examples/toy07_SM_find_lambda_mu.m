
% TESTING NEW FUNCTION cgan_l1_omega.m
% EXPERIMENT ON  LARGE OVERLAPPING BLOCKS AND SPARSE MATRIX ON OBSRVED VAR
%non oriented graph structure learning with latent variables
% 
% We build the complete model and sample from it
% We assume latept variables ndependept

%% add paths
clc
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

% %% data
run('../../toy-data/toy07.m');k=30;

%% our norm psd with decomposition S-M sparse_omega_lgm
p=po;
param.max_nb_main_loop=100;%2;%1000
param.verbose=1;
inputData.Y=-eye(po);

%%
% param.f=4;
% inputData.X1=S^.5;
% inputData.X2=S^.5;
%%
param.f=5;
param.verbose=1;
inputData.X=S;
%%
inputData.Y=-eye(po);
param.cardfun=inf*ones(1,p);
param.cardfun(k)=1;

%%% n=5000
% c=1;
% param.lambda=c*.3; %lamda ~ 2/k*mu
% param.mu=c*.1;
c=1;
param.lambda=c*.8; %lamda ~ 2/k*mu
param.mu=c*.1;

%%
%% Starting solution

ActiveSet.max_atom_count_reached=0;
ActiveSet.I={};
ActiveSet.alpha= [];
ActiveSet.atoms=pl;
ActiveSet.atom_count = pl;
if param.f==4
%     [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_sym(Doo,0);
elseif param.f==5
    [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_SM(Doo,0);
end
[ActiveSet.I_l1, ActiveSet.beta]=mat2l1index(-Doo,atoms_l1_sym);
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


startingZ.Z1=-Doo;
startingZ.Z2=Dol*Dol';



Z1=zeros(p);
nz=find(ActiveSet.beta>1e-15);
for j=nz'
    Z1=Z1+ActiveSet.beta(j)*reshape(atoms_l1_sym(:,ActiveSet.I_l1(j)),p,p);
end
Z2=zeros(p);
nz=find(ActiveSet.alpha>1e-15);
for j=nz'
    u=ActiveSet.atoms(:,j);
    Z2=Z2+ActiveSet.alpha(j)*(u*u');
end
Z=Z1+Z2;


%% blocks
[Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param,startingZ,ActiveSet);
obj_l1_om=hist.obj(end);
% nb good edges
nnz0=sum(Doo(:)~=0);
nnz=sum(Z1(:)~=0);
tp=sum(sum((Doo.*Z1<0)));
vv= (Doo==0 & Z1~=0) | (Doo.*Z1>0);
fp=sum(sum(vv));
prec=tp/(tp+fp);

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

keyboard;

%% tr+l1
% param.lambda=.3; %lamda ~ 2/k*mu
% param.mu=0.1;
param.lambda=.5; %lamda ~ 2/k*mu
param.mu=0.2;
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

%%

% 
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
% 
% 
% %% saving
% % %filename = ['lggm2_' datestr(datetime('now'),'yyyymmddTHHMMSS') '.ps'];
% % filename = ['lggm4_' datestr(clock) '.ps'];
% % %print ( '-dpsc2', filename, '-f1' )
% % print ( '-dpsc2', filename, '-append', '-f1' )
% % print ( '-dpsc2', filename, '-append', '-f2' )
% % print ( '-dpsc2', filename, '-append', '-f3' )
% % print ( '-dpsc2', filename, '-append', '-f4' )
% % print ( '-dpsc2', filename, '-append', '-f5' )
% 
% % % keyboard
% % save('lggm4_03_16','k','p','n','inputData','Dfull','Dmargo', ...
% % 'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output',...
% % 'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr');
% 
obj_l1_om
hist_tr.obj(end)

% nb good edges
nnz0=sum(Doo(:)~=0);
nnz=sum(Z1(:)~=0);
tp=sum(sum((Doo>0 & Z1<0)))+sum(sum((Doo<0 & Z1>0)));
fp=sum(sum((Doo==0 & Z1~=0)));
prec=tp/(tp+fp);