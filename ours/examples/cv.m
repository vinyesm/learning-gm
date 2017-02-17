clear all; close all; clc;

%parpool(4);

%%
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
run('../../toy-data/three_blocks_same_size.m');

%% our norm psd with decomposition S-M sparse_omega_lgm

% param.cardfun(p)=1;

%choice s.t. mu*k^=lambda
% aa=0.01;
% param.lambda=aa; %lamda ~ 2/k*mu
% param.mu=0.05;

%%
% lambda < k*mus
las=10.^linspace(0,-4,4);
pair=[];
count=1;
for i=1:length(las)
    for j=1:i
        pair(count).lambda=las(i);
        pair(count).mu=las(j);
        count=count+1;
    end
end

%%
nbfold=10;
partitions = cvpartition(n,'KFold',10);
%c=cvpartition(n,'LeaveOut');
% c.test(1) c.training(1)

parfor j=1:partitions.NumTestSets
    Xtrain=X(:,partitions.training(j));
    Xtest=X(:,partitions.test(j));
    Strain=cov(Xtrain');
    Stest=cov(Xtrain');
    for jj=1:length(pair)
        %% blocks
        [Dfin1{j}{jj}] = f1(Strain,pair(jj).lambda,pair(jj).mu,5);
        [Dfin2{j}{jj}] = f2(Strain,pair(jj).lambda,pair(jj).mu);
    end
end

figure(1);clf;
cc=1;
for j=1:partitions.NumTestSets
    for jj=1:length(pair)
        subplot(partitions.NumTestSets,length(pair),cc);
        imagesc(abs(Dfin1{j}{jj}));
        cc=cc+1;
        pbaspect([1 1 1]);
    end
end


