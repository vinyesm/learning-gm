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

objective= @(S05,Z) .5*norm(S05*Z*S05+eye(size(Z,1)),'fro')^2;

%% our norm psd with decomposition S-M sparse_omega_lgm

% param.cardfun(p)=1;

%choice s.t. mu*k^=lambda
% aa=0.01;
% param.lambda=aa; %lamda ~ 2/k*mu
% param.mu=0.05;

%%
% lambda < k*mus
jcut=3;
las=10.^linspace(0,-4,8);%0,-4,4
pair=[];
count=1;
for i=1:length(las)
    for j=max(1,i-jcut):i
        pair(count).lambda=las(i);
        pair(count).mu=las(j);
        count=count+1;
    end
end

%%
nbfold=5;
partitions = cvpartition(n,'KFold',nbfold);
%c=cvpartition(n,'LeaveOut');
% c.test(1) c.training(1)


parfor j=1:partitions.NumTestSets
% for j=1:partitions.NumTestSets
    Xtrain=X(:,partitions.training(j));
    Xtest=X(:,partitions.test(j));
    Strain=cov(Xtrain');
    Stest=cov(Xtest');
    for jj=1:length(pair)
        %% blocks
        [Dfin1{j}{jj},Z1] = f1(Strain,pair(jj).lambda,pair(jj).mu,5);
        cv1cell{j}{jj} = objective(Stest^.5,Z1); 
        [Dfin2{j}{jj},Z2] = f2(Strain,pair(jj).lambda,pair(jj).mu);
        cv2cell{j}{jj} = objective(Stest^.5,Z2);        
    end
end

cv1=zeros(partitions.NumTestSets,length(pair));
cv2=zeros(partitions.NumTestSets,length(pair));
for j=1:partitions.NumTestSets
    for jj=1:length(pair)
        cv1(j,jj)=cv1cell{j}{jj};
        cv2(j,jj)=cv2cell{j}{jj};
    end
end
%%
mcv1=mean(cv1);
mcv2=mean(cv2);
cv1grid=inf*ones(length(las));
cv2grid=inf*ones(length(las));
count=1;
mincv1=inf;
mincv2=inf;
p1=0;
p2=0;
for i=1:length(las)
    for j=max(1,i-jcut):i
        cv1grid(i,j)=mcv1(count);
        if mcv1(count)<mincv1
            mincv1=mcv1(count);
            p1=count;
        end
        cv2grid(i,j)=mcv2(count);
        if mcv2(count)<mincv2
            mincv2=mcv2(count);
            p2=count;
        end
        count=count+1;
    end
end


% figure(1);clf;
% cc=1;
% for j=1:partitions.NumTestSets
%     for jj=1:length(pair)
%         subplot(partitions.NumTestSets,length(pair),cc);
%         imagesc(abs(Dfin1{j}{jj}));
%         cc=cc+1;
%         pbaspect([1 1 1]);
%     end
% end
%%

figure(1);clf;
for j=1:partitions.NumTestSets
    subplot(2,partitions.NumTestSets,j);
    imagesc(abs(Dfin1{j}{p1}));
    pbaspect([1 1 1]);
    subplot(2,partitions.NumTestSets,partitions.NumTestSets+j);
    imagesc(abs(Dfin2{j}{p2}));
    pbaspect([1 1 1]);
end


%% chosen pair,

% p1=24
% lambda: 3.7276e-04
% mu: 0.0720

%p1=24;
k=5;
lambda=pair(p1).lambda;
mu=pair(p1).mu;

S=cov(X');
p=size(S,1);
param.f=4;
param.diag=0;
param.PSD=true;
param.max_nb_main_loop=100;%2;%1000
param.powerIter=500;
param.stPtPowerIter=1000;
param.niterPS=10000;%5000
param.epsStop=1e-8;
param.PSdualityEpsilon=1e-8;
param.k=0;
param.PSmu=0; %strong convexity
param.verbose=1;
param.debug=0;
param.sloppy=0;
param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
inputData.Y=-eye(p);
inputData.X1=S^.5;
inputData.X2=S^.5;
param.cardfun=inf*ones(1,p);
param.cardfun(k)=1;
param.lambda=lambda;
param.mu=mu;
[Zf Z1f Z2f ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
if ~isempty(ActiveSet.alpha)
    Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
    nl=size(ActiveSet.atoms,2);
    Dfin=zeros(p+nl);
    Dfin(1:nl,1:nl)=eye(nl);
    Dfin((nl+1):(nl+p),(nl+1):(nl+p))=-Z1f;
    Dfin(1:nl,(nl+1):(nl+p))=Uso';
    Dfin((nl+1):(nl+p),1:nl)=Uso;
else
    Dfin=Z1f;
end

figure(2);clf;
imagesc(abs(Dfin));
pbaspect([1 1 1]);





