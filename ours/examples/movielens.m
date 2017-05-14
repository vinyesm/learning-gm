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

run('pp_movielens.m');
% run('pp_movielens_2.m');
k=20;
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


%% Starting solution

% Doo=eye(60);
% ActiveSet = {};
% if param.f==4
%     [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_sym(Doo,0);
% elseif param.f==5
%     [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_SM(Doo,0);
% end
% [ActiveSet.I_l1, ActiveSet.beta]=mat2l1index(-Doo,atoms_l1_sym);
% 
% ActiveSet.I = {};
% ActiveSet.U = {};
% ActiveSet.Sigma = {};
% ActiveSet.Z = {};
% ActiveSet.tracenorm = {};
% ActiveSet.fronorm = {};
% ActiveSet.k = {};
% ActiveSet.atomsSupport = {};
% ActiveSet.alpha= [];
% ActiveSet.atom_count = 0;
% ActiveSet.atoms=[];
% ActiveSet.max_atom_count_reached=0;
% 
% startingZ.Z1=-Doo;
% startingZ.Z2=zeros(60);
% 
% Z1=zeros(p);
% nz=find(ActiveSet.beta>1e-15);
% for j=nz'
%     Z1=Z1+ActiveSet.beta(j)*reshape(atoms_l1_sym(:,ActiveSet.I_l1(j)),p,p);
% end
% Z2=zeros(p);
% nz=find(ActiveSet.alpha>1e-15);
% for j=nz'
%     u=ActiveSet.atoms(:,j);
%     Z2=Z2+ActiveSet.alpha(j)*(u*u');
% end
% Z=Z1+Z2;

%%
%% Starting solution
param.cardfun=inf*ones(1,p);
param.cardfun(k)=1;

Doo=eye(60);
pl=3;
ks=[20 20 20];
Dol=zeros(60,pl);
Dol(1:20,1)=ones(20,1)/sqrt(20);
Dol(21:40,2)=ones(20,1)/sqrt(20);
Dol(41:60,3)=ones(20,1)/sqrt(20);

ActiveSet.max_atom_count_reached=0;
ActiveSet.I={};
ActiveSet.alpha= [];
ActiveSet.atoms=pl;
ActiveSet.atom_count = pl;
if param.f==4
     [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_sym(Doo,0);
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

%% reg param
% lam .4 mu .1 f==4

lam=.5;
gam=.1;
% lam=.5;
% gam=.1;
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
% [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
[Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param,startingZ, ActiveSet);
% obj_l1_om=hist.obj(end);

save(['ml_' num2str(param.f) '_' num2str(param.mu) '_' num2str(param.lambda) ],'param', ...
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

figure(2);clf;
subplot(1,4,1)
imagesc(abs(Dfin));
pbaspect([1 1 1]);
subplot(1,4,2)
imagesc(abs(Dfin)>0);
pbaspect([1 1 1]);
subplot(1,4,3)
imagesc(abs(Z));
pbaspect([1 1 1]);
subplot(1,4,4)
imagesc(abs(Z2));
pbaspect([1 1 1]);

keyboard;

%% tr+l1
%% Starting solution
param.cardfun=inf*ones(1,p);
param.cardfun(p)=1;

% Doo=eye(60);
% pl=3;
% ks=[60 60 60];
% Dol=zeros(60,pl);
% Dol(1:20,1)=ones(20,1)/sqrt(20);
% Dol(21:40,2)=ones(20,1)/sqrt(20);
% Dol(41:60,3)=ones(20,1)/sqrt(20);

Doo=-Z1;
Dol=Uso;
pl=size(Uso,2);
ks=60*ones(1,size(Uso,2));

ActiveSet.max_atom_count_reached=0;
ActiveSet.I={};
ActiveSet.alpha= [];
ActiveSet.atoms=pl;
ActiveSet.atom_count = pl;
if param.f==4
     [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_sym(Doo,0);
elseif param.f==5
    [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_SM(Doo,0);
end
[ActiveSet.I_l1, ActiveSet.beta]=mat2l1index(-Doo,atoms_l1_sym);
ActiveSet.k=mat2cell(ks,1,ones(1,length(ks)));
ActiveSet.alpha=sum(Dol.^2)';
ActiveSet.atoms=sparse(bsxfun(@rdivide, Dol, sqrt(sum(Dol.^2))));
for i=1:pl
    ActiveSet.I{i}=(1:p)';
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


param.lambda=.6; %lamda ~ 2/k*mu
param.mu=.1;
param.max_nb_main_loop=2;
param.niterPS=10000;
% param.cardfun=inf*ones(1,p);
% param.cardfun(p)=1;
[Z_tr Z1_tr Z2_tr ActiveSet_tr hist_tr param_tr flag_tr output_tr] = cgan_l1_omega(inputData,param);
% [Z_tr Z1_tr Z2_tr ActiveSet_tr hist_tr param_tr flag_tr output_tr] = cgan_l1_omega(inputData,param,startingZ,ActiveSet);

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

save(['ml_' num2str(param.f) '_' num2str(param.mu) '_' num2str(param.lambda) ],'param', ...
'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr');
% save('ml','k','p','n','inputData','Dfull','Dmargo', ...
% 'Z', 'Z1', 'Z2', 'ActiveSet', 'hist' ,'param', 'flag' ,'output',...
% 'Z_tr', 'Z1_tr', 'Z2_tr', 'ActiveSet_tr', 'hist_tr', 'param_tr', 'flag_tr', 'output_tr');

% obj_l1_om
% hist_tr.obj(end)