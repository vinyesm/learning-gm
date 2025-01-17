% Toy example
% 3 blocks of size 20x10x5
% S is a chain 

clear all; clc;
addpath('../ours/utils');

ks=[20 15 10];
pl=length(ks);
po=sum(ks);
n=100*po;

c=0;
pt=pl+po;
clo=.7./sqrt(ks);
% clo=.5./(ones(1,pl));
corr=zeros(po);

corr_all= zeros(pl+po,pl+po);
corr_all(1:pl,1:pl)=eye(pl);

for i=1:po
    corr(i,i)=1;
    if i<po
        corr(i,i+1)=c;
        corr(i+1,i)=c;
    end
end
corr_all(pl+1:end,pl+1:end)=corr;

for i=1:pl
    if i>1
        i1=sum(ks(1:(i-1)))+1+pl;
    else
        i1=1+pl;
    end
    i2=i1+ks(i)-1;
%     Uall((pl+i1):(pl+i2),i)=-.5*ones(ks(i),1);
    corr_all(i1:i2,i)=ones(ks(i),1)*clo(i);
    corr_all(i,i1:i2)=corr_all(i1:i2,i)';
end

%%
%% observed correlation matrix
CostMatrix = rand(po);
[ Tree,Cost ] =  UndirectedMaximumSpanningTree ( CostMatrix );
Tree=.2*Tree;
Tree = Tree - diag(diag(Tree)) + eye(po);
corr_all(pl+1:end,pl+1:end)=Tree;


figure(1); clf;
imagesc(corr_all);
pbaspect([1 1 1]);

%%
Dfull=corr_all;
Doo=corr_all(pl+1:end,pl+1:end);
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

%% active set
ActiveSet.max_atom_count_reached=0;
ActiveSet.I={};
ActiveSet.alpha= [];
ActiveSet.atoms=pl;
ActiveSet.atom_count = pl;
[ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_sym(Doo,0);
[ActiveSet.I_l1, ActiveSet.beta]=mat2l1index(-Doo,atoms_l1_sym);
ActiveSet.k=mat2cell(ks,1,ones(1,length(ks)));
ActiveSet.alpha=sum(Dol.^2)';
ActiveSet.atoms=sparse(bsxfun(@rdivide, Dol, sqrt(sum(Dol.^2))));
for i=1:pl
    ActiveSet.I{i}=find(Dol(:,i));
end


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

