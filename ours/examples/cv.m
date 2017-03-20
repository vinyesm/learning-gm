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
% run('../../toy-data/three_blocks_same_size.m');k=5;rank=3;p=size(X,1);Z0=eye(p);
run('../../toy-data/three_large_blocks_same_size.m');k=10; rank=5;p=size(X,1);Z0=eye(p);

objective = @(S05,Z) .5*norm(S05*Z*S05+eye(size(Z,1)),'fro')^2;
rankdiff = @(ActiveSet) abs((rank-size(ActiveSet.atoms,2)));
suppdiff = @(Z,Z0) (sum(sign(Z(:))==sign(Z0(:))));

%% our norm psd with decomposition S-M sparse_omega_lgm

% param.cardfun(p)=1;

%%
% lambda < k*mus
% jcut=inf;
mus=[0.1];
las=[0.1 0.5 0.01 0.005 0.001];
%las=10.^linspace(0,-4,8);%0,-4,4
pair=[];
count=1;
for i=1:length(las)
    for j=1:length(mus)
        if mus(j)>=las(i),
            pair(count).lambda=las(i);
            pair(count).mu=mus(j);
            count=count+1;
        end
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
        [Dfin1{j}{jj},Z1,Z11,ActiveSet] = f1(Strain,pair(jj).lambda,pair(jj).mu,k);
        cv1cell{j}{jj} = objective(Stest^.5,Z1); 
        cv1cell_supp{j}{jj} = suppdiff(-Z11,Z0);
        cv1cell_rank{j}{jj} = rankdiff(ActiveSet);
    end
end

save('cv01midf1', 'k','X', 'pair', 'partitions', 'cv1cell' ,'Dfin1');
%%
parfor j=1:partitions.NumTestSets
% for j=1:partitions.NumTestSets
    Xtrain=X(:,partitions.training(j));
    Xtest=X(:,partitions.test(j));
    Strain=cov(Xtrain');
    Stest=cov(Xtest');
    for jj=1:length(pair)
        [Dfin2{j}{jj},Z2,Z21,ActiveSet] = f2(Strain,pair(jj).lambda,pair(jj).mu);
        cv2cell{j}{jj} = objective(Stest^.5,Z2);
        cv2cell_supp{j}{jj} = suppdiff(-Z21,Z0);
        cv2cell_rank{j}{jj} = rankdiff(ActiveSet);
    end
end

save('cv01midf2', 'k','X',  'pair', 'partitions','Dfin1', 'Dfin2', 'cv1cell' ,'cv2cell');

[cv1_obj, p1_obj, mincv1_obj]= cv_out(cv1cell,partitions.NumTestSets,pair,las,mus);
[cv2_obj, p2_obj, mincv2_obj]= cv_out(cv2cell,partitions.NumTestSets,pair,las,mus);
[cv1_supp, p1_supp, mincv1_supp]= cv_out(cv1cell_supp,partitions.NumTestSets,pair,las,mus);
[cv2_supp, p2_supp, mincv2_supp]= cv_out(cv2cell_supp,partitions.NumTestSets,pair,las,mus);
[cv1_rank, p1_rank, mincv1_rank]= cv_out(cv1cell_rank,partitions.NumTestSets,pair,las,mus);
[cv2_rank, p2_rank, mincv2_rank]= cv_out(cv2cell_rank,partitions.NumTestSets,pair,las,mus);

%%
figure(5); clf;
subplot(1,2,1)
plot(cv1_obj(:,1)./norm(cv1_obj(:,1)),'k'); hold on;
plot(cv1_supp(:,1)./norm(cv1_supp(:,1)),'r'); hold on;
plot(cv1_rank(:,1)./norm(cv1_rank(:,1)),'b');
legend('obj','supp','rank');
title(['lambda selection l1+om for mu=' num2str(mus(1)) ]);
subplot(1,2,2)
plot(cv2_obj(:,1)./norm(cv2_obj(:,1)),'k'); hold on;
plot(cv2_supp(:,1)./norm(cv2_supp(:,1)),'r'); hold on;
plot(cv1_rank(:,1)./norm(cv2_rank(:,1)),'b');
legend('obj','supp','rank');
title(['lambda selection l1+tr for mu=' num2str(mus(1)) ]);
% legend('show','Location','southwest');

%%
% p1=p1_obj;
% p2=p2_obj;
p1=p1_rank;
p2=p2_rank;
% p1=p1_supp;
% p2=p2_supp;

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

% p1=p1_obj;
% p2=p2_obj;
p1=p1_rank;
p2=p2_rank;
% p1=p1_supp;
% p2=p2_supp;



%%
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
param.niterPS=500;%5000
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
    Dfin_l1_om=zeros(p+nl);
    Dfin_l1_om(1:nl,1:nl)=eye(nl);
    Dfin_l1_om((nl+1):(nl+p),(nl+1):(nl+p))=-Z1f;
    Dfin_l1_om(1:nl,(nl+1):(nl+p))=Uso';
    Dfin_l1_om((nl+1):(nl+p),1:nl)=Uso;
else
    Dfin_l1_om=Z1f;
end

figure(2);clf;
imagesc(abs(Dfin_l1_om));
pbaspect([1 1 1]);

%%
param.max_nb_main_loop=2;%2;%1000
lambda=pair(p2).lambda;
mu=pair(p2).mu;
param.lambda=lambda;
param.mu=mu;
param.cardfun=inf*ones(1,p);
param.cardfun(p)=1;
[Zff Z1ff Z2ff ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
if ~isempty(ActiveSet.alpha)
    Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
    nl=size(ActiveSet.atoms,2);
    Dfin_l1_tr=zeros(p+nl);
    Dfin_l1_tr(1:nl,1:nl)=eye(nl);
    Dfin_l1_tr((nl+1):(nl+p),(nl+1):(nl+p))=-Z1ff;
    Dfin_l1_tr(1:nl,(nl+1):(nl+p))=Uso';
    Dfin_l1_tr((nl+1):(nl+p),1:nl)=Uso;
else
    Dfin_l1_tr=Z1f;
end


figure(3);clf;
imagesc(abs(Dfin_l1_tr));
pbaspect([1 1 1]);


%%
% save('cv01', 'p', 'k', 'inputData','X', 'pair', ...
%     'cv1_obj', 'p1_obj', 'mincv1_obj',...
%     'cv2_obj', 'p2_obj', 'mincv2_obj',...
%     'cv1_supp', 'p1_supp', 'mincv1_supp',...
%     'cv2_supp', 'p2_supp', 'mincv2_supp',...
%     'cv1_rank', 'p1_rank', 'mincv1_rank',...
%     'cv2_rank', 'p2_rank', 'mincv2_rank',...
%     'Dfin1', 'Dfin2');
save('cv02', 'p', 'k', 'inputData','X', 'pair', ...
    'cv1_obj', 'p1_obj', 'mincv1_obj',...
    'cv2_obj', 'p2_obj', 'mincv2_obj',...
    'cv1_supp', 'p1_supp', 'mincv1_supp',...
    'cv2_supp', 'p2_supp', 'mincv2_supp',...
    'cv1_rank', 'p1_rank', 'mincv1_rank',...
    'cv2_rank', 'p2_rank', 'mincv2_rank',...
    'Dfin1', 'Dfin2');
