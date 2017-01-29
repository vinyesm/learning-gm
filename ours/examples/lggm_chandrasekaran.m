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
n=1000; % number of samples

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


descr = {'Plot of the function:';
    'y = A{\ite}^{-\alpha{\itt}}';
    ' ';
    'With the values:';
    'A = 0.25';
    '\alpha = .005';
    't = 0:1000'};

descr = {'Plot of the groundtruth concentration matrix';
    'of the complete model'};

figure(1);clf;
subplot(1,2,1)
axis off;
text(0.5,.5,descr)
pbaspect([1 1 1]);
subplot(1,2,2)
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
% S=inv(Dmargo);

%% plotting

figure(2);clf;
subplot(2,2,1);
imagesc(abs(Dmargo));
pbaspect([1 1 1]);
title('true marginal conc. mat.');
subplot(2,2,2);
imagesc(abs(inv(S)));
pbaspect([1 1 1]);
title('observed conc. mat.');
subplot(2,2,3);
imagesc(abs(inv(Dmargo)));
pbaspect([1 1 1]);
title('true cov');
subplot(2,2,4);
imagesc(abs(S));
pbaspect([1 1 1]);
title('observed cov');


%% LVGGM Chandrasekaran S-L, (Sparse-Low Rank)


HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

addpath LogdetPPA-0/;
addpath LogdetPPA-0/solver/;
addpath LogdetPPA-0/solver/mexfun/;
addpath LogdetPPA-0/util/;
ttime  = clock;