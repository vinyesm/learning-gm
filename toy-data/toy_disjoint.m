%% TOY EXAMPLE WITH TREE STRUCTURE ON OBSERVED VARIABLES
clear all; clc;
addpath ../ours/utils/

%%

ks=[15 15 15];
pl=length(ks);
po=sum(ks);
n=50*po;

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
Tree=.2*Tree;
Tree = Tree - diag(diag(Tree)) + eye(po);
corr_all(pl+1:end,pl+1:end)=Tree;

%% observed-latent
for i=1:pl
    if i>1
        i1=sum(ks(1:(i-1)))+1+pl;
    else
        i1=1+pl;
    end
    i2=i1+ks(i)-1;
    corr_all(i1:i2,i)=ones(ks(i),1)*clo(i);
    corr_all(i,i1:i2)=ones(ks(i),1)'*clo(i);
end


figure(1); clf;
imagesc(corr_all);
pbaspect([1 1 1]);

%%
Dfull=corr_all;
Doo=Tree;
Dol=corr_all(pl+1:end,1:pl);
Dmargo=Doo-Dol*Dol';


%%
descr = {'Plot of the groundtruth '
    'concentration matrix';
    'of the complete model'};

figure(1);clf;
subplot(1,2,1)
axis off;
text(0,.5,descr)
pbaspect([1 1 1]);
subplot(1,2,2)
imagesc(abs(Dfull));
colorbar;
pbaspect([1 1 1]);
title('true complete conc. mat.');

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


 
