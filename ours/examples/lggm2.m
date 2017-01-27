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
addpath('../other');
addpath('../prox');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

%% settings
clear all; clc; close all;
n=5000; % number of samples

pl=3; % number of latept variables
po=15;% number of observed variables
pt=po+pl;

slo=.2; % level of sparsity on latent-observed block of conceptration mat
soo=.2; % level of sparsity on latent-observed block of conceptration mat
Too=0*tril(rand(po,po)<soo,-1);
Doo=.2*(Too+Too')+eye(po);
Dol=zeros(po,pl);
Dol(1:5,1)=ones(5,1);
Dol(6:10,2)=ones(5,1);
Dol(11:15,3)=ones(5,1);
Dol=.2*Dol;

%% construction of the conceptration matrix
Dfull=zeros(pt);
Dfull(1:pl,1:pl)=eye(pl);
Dfull((pl+1):pt,(pl+1):pt)=Doo;
Dfull(1:pl,(pl+1):pt)=Dol';
Dfull((pl+1):pt,1:pl)=Dol;

Dfull=.5*(Dfull+Dfull');
Dmargo=Doo-Dol*Dol';



figure(1);clf;
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');


if eigs(Dfull,1,'sa')<0
    error('The conceptration matrix of the complete model is not PSD \n')
end

%% sampling data
mu=zeros(1,pt); % vector of means
Xfull=mvnrnd(mu, inv(Dfull), n)';
X=Xfull((pl+1):pt,:);
S=cov(X');
S=inv(Dmargo);

%% plotting

figure(2);clf;
subplot(1,2,1);
imagesc(abs(Dmargo));
pbaspect([1 1 1]);
title('true marginal conc. mat.');
subplot(1,2,2);
imagesc(abs(inv(S)));
pbaspect([1 1 1]);
title('observed conc. mat.');

%% our norm psd with decomposition S-M
%% sparse_omega_lgm;
p=po;
inputData.X1=S^.5;
param.mu=1e-1;
param.lambda=1e-2;
param.rho=1;
param.cardfun=inf*ones(1,p);
param.cardfun(5)=1;
% param.cardfun=(1:(p)).^.2;
[Aso,Mso,Sso,Eso,Mso_as,out] = sparse_omega_lgm( inputData, param);

obj1=out.obj;
loss1=.5*norm(S^.5*(Sso-Mso)*S^.5-eye(po),'fro')^2;
pens1=param.mu*sum(abs(Sso(:))>0);
peno1=param.lambda*out.Mnorm;

obj0=.5*norm(S^.5*Dmargo*S^.5-eye(po),'fro')^2+param.mu*sum(abs(Doo(:))>0)+param.lambda*pl*5*0.2;
loss0=.5*norm(S^.5*Dmargo*S^.5-eye(po),'fro')^2;
pens0=param.mu*sum(abs(Doo(:))>0);
peno0=param.lambda*pl*5*0.2;

fprintf('ground truth sol : obj=%f  loss=%f  pens=%f   peno=%f\n',obj0,loss0,pens0,peno0);
fprintf('our  sol         : obj=%f  loss=%f  pens=%f   peno=%f\n',obj1,loss1,pens1,peno1);

if norm(Mso,'fro')>1e-16
    Uso=bsxfun(@times,sqrt(Mso_as.alpha)',Mso_as.atoms);
    nl=length(Mso_as.alpha);
    Dso=zeros(p+nl);
    Dso(1:nl,1:nl)=eye(nl);
    Dso((nl+1):(nl+p),(nl+1):(nl+p))=Sso;
    Dso(1:nl,(nl+1):(nl+p))=Uso';
    Dso((nl+1):(nl+p),1:nl)=Uso;
else
    Dso=Sso;
end

figure(3);clf;
subplot(2,2,1)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
subplot(2,2,2)
imagesc(abs(Dso));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
subplot(2,2,3)
imagesc(abs(Dfull)>1e-15);
pbaspect([1 1 1]);
title('true support');
subplot(2,2,4)
imagesc(abs(Dso)>1e-15);
pbaspect([1 1 1]);
title('estimated support');

figure(4);clf;
subplot(2,2,1);
imagesc(abs(Dmargo));
pbaspect([1 1 1]);
title('true marginal conc. mat.');
subplot(2,2,2);
imagesc(abs(Sso-Mso));
pbaspect([1 1 1]);
title('observed conc. mat.');
subplot(2,2,3)
imagesc(abs(Dmargo)>1e-15);
pbaspect([1 1 1]);
title('true support');
subplot(2,2,4)
imagesc(abs(Sso-Mso)>1e-15);
pbaspect([1 1 1]);
title('estimated support');



