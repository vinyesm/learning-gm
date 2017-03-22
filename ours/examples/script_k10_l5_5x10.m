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
mus=[0.1];
las=[0.1 0.5 0.01 0.005 0.001];
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

parfor j=1:partitions.NumTestSets
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

save('cv_large_blocks_midf1', 'k','X', 'pair', 'partitions', 'cv1cell' ,'Dfin1');
%%
parfor j=1:partitions.NumTestSets
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

save('cv_large_blocks_midf2', 'k','X',  'pair', 'partitions','Dfin1', 'Dfin2', 'cv1cell' ,'cv2cell');

[cv1_obj, p1_obj, mincv1_obj]= cv_out(cv1cell,partitions.NumTestSets,pair,las,mus);
[cv2_obj, p2_obj, mincv2_obj]= cv_out(cv2cell,partitions.NumTestSets,pair,las,mus);
[cv1_supp, p1_supp, mincv1_supp]= cv_out(cv1cell_supp,partitions.NumTestSets,pair,las,mus);
[cv2_supp, p2_supp, mincv2_supp]= cv_out(cv2cell_supp,partitions.NumTestSets,pair,las,mus);
[cv1_rank, p1_rank, mincv1_rank]= cv_out(cv1cell_rank,partitions.NumTestSets,pair,las,mus);
[cv2_rank, p2_rank, mincv2_rank]= cv_out(cv2cell_rank,partitions.NumTestSets,pair,las,mus);


save('cv_large_blocks', 'p', 'k','X', 'pair', ...
    'cv1_obj', 'p1_obj', 'mincv1_obj',...
    'cv2_obj', 'p2_obj', 'mincv2_obj',...
    'cv1_supp', 'p1_supp', 'mincv1_supp',...
    'cv2_supp', 'p2_supp', 'mincv2_supp',...
    'cv1_rank', 'p1_rank', 'mincv1_rank',...
    'cv2_rank', 'p2_rank', 'mincv2_rank',...
    'Dfin1', 'Dfin2',...
    'las','mus');
