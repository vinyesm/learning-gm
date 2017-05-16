%% TOY EXAMPLE WITH TREE STRUCTURE ON OBSERVED VARIABLES
%% AND OVERLAPPING BLOCKS
%%
clear all; clc;
addpath ../ours/utils/

%%
lsupp=[1:30;21:50;41:70;61:90;81:110;101:130;121:150];
pl=size(lsupp,1);
ks=30*ones(1,pl);
po=max(lsupp(:));
%n=100*po;
n=10*po;

pt=pl+po;
clo=.7./sqrt(ks);
% clo=.5./(ones(1,pl));
corr=zeros(po);

%% correlation matrix
corr_all= zeros(pl+po,pl+po);
corr_all(1:pl,1:pl)=eye(pl);

%% observed correlation matrix
CostMatrix = rand(po);
[ Tree,Cost ] =  UndirectedMaximumSpanningTree ( CostMatrix );
Tree=Tree;
Tree = Tree - diag(diag(Tree)) + eye(po);
corr_all(pl+1:end,pl+1:end)=Tree;

%% observed-latent
for i=1:pl
    corr_all(pl+lsupp(i,:),i)=ones(ks(i),1);
    corr_all(i,pl+lsupp(i,:))=ones(ks(i),1)';
end


figure(1); clf;
imagesc(corr_all);
pbaspect([1 1 1]);

%%
eta=1.05;
G=corr_all;
L=eta*diag(sum(G))-G;
LAM=diag(diag(inv(L)));
LAM=LAM.^.5;
Dfull=LAM*L*LAM;
Dfull=.5*(Dfull+Dfull');
Doo=Dfull(pl+1:end,pl+1:end);
Dol=Dfull(pl+1:end,1:pl);
Dmargo=Doo-Dol*(inv(Dfull(1:pl,1:pl)))*Dol';
Dmargo=.5*(Dmargo+Dmargo');

%%

figure(1);clf;
subplot(1,2,1)
imagesc(abs(Dfull));
axis square
subplot(1,2,2)
imagesc(abs(Dmargo));
axis square

if eigs(Dfull,1,'sa')<0
    error('The concentration matrix of the complete model is not PSD \n')
end
if eigs(Dmargo,1,'sa')<0
    error('Dmargo matrix of the complete model is not PSD \n')
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


 
