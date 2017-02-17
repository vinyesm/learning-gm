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
run('../../toy-data/three_large_blocks_same_size.m');

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

parfor j=1:length(partitions)
    %%
    p=po;
    param.f=4;
    param.diag=0;
    param.PSD=true;
    param.max_nb_main_loop=2;%100;%2;%1000
    param.powerIter=500;
    param.stPtPowerIter=1000;
    param.niterPS=5;%5000;%10000;%5000
    param.epsStop=1e-8;
    param.PSdualityEpsilon=1e-8;
    param.k=0;
    param.PSmu=0; %strong convexity
    param.verbose=1;
    param.debug=0;
    param.sloppy=0;
    param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
    inputData.Y=-eye(po);
    param.cardfun=inf*ones(1,p);
    param.cardfun(k)=1;
    %%
    Xtrain=X(:,part.training(j));
    Xtest=X(:,part.test(j));
    Strain=cov(Xtrain');
    Stest=cov(Xtrain');
    inputData.X1=Strain^.5;
    inputData.X2=Strain^.5;
    for jj=1:length(pairs)
        param.lambda=pairs(jj).lambda;
        param.mu=pairs(jj).mu;
        %% blocks
        [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
    end
end