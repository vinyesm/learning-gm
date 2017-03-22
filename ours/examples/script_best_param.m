%% SEARCHING BEST PAIR lambda, mu (knowing ground truth)
%%
clear all; close all; clc;

delete(gcp)
parpool(4);

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
% run('../../toy-data/three_blocks_same_size.m');k=5;rank=3;p=size(X,1);Z0=eye(p);S=cov(X');
run('../../toy-data/three_large_blocks_same_size.m');k=10; rank=5;p=size(X,1);Z0=eye(p);S=cov(X');


%%
objective = @(S05,Z) .5*norm(S05*Z*S05+eye(size(Z,1)),'fro')^2;
rankdiff = @(ActiveSet) abs((rank-size(ActiveSet.atoms,2)));
suppdiff = @(Z,Z0) (sum(sign(Z(:))==sign(Z0(:))));

%%
% lambda < k*mus
c=sqrt(p/n);
mus= linspace(.01,1,10);
las=linspace(.01,1,10);
pair=[];
count=1;
for i=1:length(las)
    for j=1:length(mus)
%         if mus(j)>=las(i),
            pair(count).lambda=las(i);
            pair(count).mu=mus(j);
            count=count+1;
%         end
    end
end

%%
cv1_obj=zeros(length(pair),1);
cv1_supp=zeros(length(pair),1);
cv1_rank=zeros(length(pair),1);
cv1_lat=inf*ones(length(pair),1);
cv2_obj=zeros(length(pair),1);
cv2_supp=zeros(length(pair),1);
cv2_rank=zeros(length(pair),1);
cv2_lat=inf*ones(length(pair),1);


parfor jj=1:length(pair)
    [Dfin1{jj},Z1{jj},Z11,ActiveSet1,Uso1] = f1(S,pair(jj).lambda,pair(jj).mu,k);
    cv1_obj(jj)  = objective(S^.5,Z1{jj}); 
    cv1_supp(jj) = suppdiff(-Z11,Z0)/p^2;
    cv1_rank(jj) = rankdiff(ActiveSet1);
    if cv1_rank(jj)==0
        cv1_lat(jj) =  suppdiff(Uso1,Dol);
    end
end
parfor jj=1:length(pair)
    [Dfin2{jj},Z2{jj},Z21,ActiveSet2,Uso2] = f2(S,pair(jj).lambda,pair(jj).mu);
    cv2_obj(jj)  = objective(S^.5,Z2{jj}); 
    cv2_supp(jj) = suppdiff(-Z21,Z0)/p^2;
    cv2_rank(jj) = rankdiff(ActiveSet2);
    if cv2_rank(jj)==0
        cv2_lat(jj) =  suppdiff(Uso2,Dol);
    end
end

%%
[vo1 po1]=min(cv1_obj);
[vs1 ps1]=min(cv1_supp);
[vr1 pr1]=min(cv1_rank);
[vl1 pl1]=min(cv1_lat);

[vo2 po2]=min(cv2_obj);
[vs2 ps2]=min(cv2_supp);
[vr2 pr2]=min(cv2_rank);
[vl2 pl2]=min(cv2_lat);

%% plot
eps=1e-12;

figure(1);clf
subplot(3,2,1);
imagesc(abs(Z1{po1})>eps);
pbaspect([1 1 1])
title(['obs l1+om  best obj lam=' num2str(pair(po1).lambda) ' mu=' num2str(pair(po1).mu)] )
subplot(3,2,2);
imagesc(abs(Z2{po2})>eps);
pbaspect([1 1 1])
title(['obs l1+tr best obj  lam=' num2str(pair(po2).lambda) ' mu=' num2str(pair(po2).mu)] )
subplot(3,2,3);
imagesc(abs(Z1{ps1})>eps);
pbaspect([1 1 1])
title(['obs l1+om  best supp lam=' num2str(pair(ps1).lambda) ' mu=' num2str(pair(ps1).mu)] )
subplot(3,2,4);
imagesc(abs(Z2{ps2})>eps);
pbaspect([1 1 1])
title(['obs l1+tr best supp  lam=' num2str(pair(ps2).lambda) ' mu=' num2str(pair(ps2).mu)] )
subplot(3,2,5);
imagesc(abs(Z1{pr1})>eps);
pbaspect([1 1 1])
title(['obs l1+om  best rank lam=' num2str(pair(pr1).lambda) ' mu=' num2str(pair(pr1).mu)] )
subplot(3,2,6);
imagesc(abs(Z2{pr2})>eps);
pbaspect([1 1 1])
title(['obs l1+tr best rank  lam=' num2str(pair(pr2).lambda) ' mu=' num2str(pair(pr2).mu)] )


figure(2);clf
subplot(3,2,1);
imagesc(abs(Dfin1{po1})>eps);
pbaspect([1 1 1])
title(['lat+obs l1+om  best obj lam=' num2str(pair(po1).lambda) ' mu=' num2str(pair(po1).mu)] )
subplot(3,2,2);
imagesc(abs(Dfin2{po2})>eps);
pbaspect([1 1 1])
title(['lat+obs l1+tr best obj  lam=' num2str(pair(po2).lambda) ' mu=' num2str(pair(po2).mu)] )
subplot(3,2,3);
imagesc(abs(Dfin1{ps1})>eps);
pbaspect([1 1 1])
title(['lat+obs l1+om  best supp lam=' num2str(pair(ps1).lambda) ' mu=' num2str(pair(ps1).mu)] )
subplot(3,2,4);
imagesc(abs(Dfin2{ps2})>eps);
pbaspect([1 1 1])
title(['obs l1+tr best supp  lam=' num2str(pair(ps2).lambda) ' mu=' num2str(pair(ps2).mu)] )
subplot(3,2,5);
imagesc(abs(Dfin1{pr1})>eps);
pbaspect([1 1 1])
title(['lat+obs l1+om  best rank lam=' num2str(pair(pr1).lambda) ' mu=' num2str(pair(pr1).mu)] )
subplot(3,2,6);
imagesc(abs(Dfin2{pr2})>eps);
pbaspect([1 1 1])
title(['lat+obs l1+tr best rank  lam=' num2str(pair(pr2).lambda) ' mu=' num2str(pair(pr2).mu)] )

figure(3);clf
subplot(3,2,1);
imagesc(abs(Z1{po1}));
pbaspect([1 1 1])
title(['obs l1+om  best obj lam=' num2str(pair(po1).lambda) ' mu=' num2str(pair(po1).mu)] )
subplot(3,2,2);
imagesc(abs(Z2{po2}));
pbaspect([1 1 1])
title(['obs l1+tr best obj  lam=' num2str(pair(po2).lambda) ' mu=' num2str(pair(po2).mu)] )
subplot(3,2,3);
imagesc(abs(Z1{ps1}));
pbaspect([1 1 1])
title(['obs l1+om  best supp lam=' num2str(pair(ps1).lambda) ' mu=' num2str(pair(ps1).mu)] )
subplot(3,2,4);
imagesc(abs(Z2{ps2}));
pbaspect([1 1 1])
title(['obs l1+tr best supp  lam=' num2str(pair(ps2).lambda) ' mu=' num2str(pair(ps2).mu)] )
subplot(3,2,5);
imagesc(abs(Z1{pr1}));
pbaspect([1 1 1])
title(['obs l1+om  best rank lam=' num2str(pair(pr1).lambda) ' mu=' num2str(pair(pr1).mu)] )
subplot(3,2,6);
imagesc(abs(Z2{pr2}));
pbaspect([1 1 1])
title(['obs l1+tr best rank  lam=' num2str(pair(pr2).lambda) ' mu=' num2str(pair(pr2).mu)] )


figure(4);clf
subplot(3,2,1);
imagesc(abs(Dfin1{po1}));
pbaspect([1 1 1])
title(['lat+obs l1+om  best obj lam=' num2str(pair(po1).lambda) ' mu=' num2str(pair(po1).mu)] )
subplot(3,2,2);
imagesc(abs(Dfin2{po2}));
pbaspect([1 1 1])
title(['lat+obs l1+tr best obj  lam=' num2str(pair(po2).lambda) ' mu=' num2str(pair(po2).mu)] )
subplot(3,2,3);
imagesc(abs(Dfin1{ps1}));
pbaspect([1 1 1])
title(['lat+obs l1+om  best supp lam=' num2str(pair(ps1).lambda) ' mu=' num2str(pair(ps1).mu)] )
subplot(3,2,4);
imagesc(abs(Dfin2{ps2}));
pbaspect([1 1 1])
title(['obs l1+tr best supp  lam=' num2str(pair(ps2).lambda) ' mu=' num2str(pair(ps2).mu)] )
subplot(3,2,5);
imagesc(abs(Dfin1{pr1}));
pbaspect([1 1 1])
title(['lat+obs l1+om  best rank lam=' num2str(pair(pr1).lambda) ' mu=' num2str(pair(pr1).mu)] )
subplot(3,2,6);
imagesc(abs(Dfin2{pr2}));
pbaspect([1 1 1])
title(['lat+obs l1+tr best rank  lam=' num2str(pair(pr2).lambda) ' mu=' num2str(pair(pr2).mu)] )


