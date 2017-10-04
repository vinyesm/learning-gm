%% from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13204

addpath ../../MILE/

clc; clear all;
data=importdata('GSE13204-GPL7473_series_matrix.txt');
td=data.textdata;
% data=importdata('GSE13204_family.soft');
% data=tdfread('GSE13204-GPL7473_series_matrix.txt');

%%
genes=data.textdata(66:1522,1);
ngen=length(genes);

% important genes
% https://genecards.weizmann.ac.il/cgi-bin/geneannot/GA_search
imp_genes={ '206674_at', '221691_x_at', '221923_s_at','204039_at','205051_s_at','202647_s_at',...
    '1559856_s_at','212079_s_at','212078_s_at','212076_at', '1565436_s_at'...
    '206067_s_at'};
imp_genes_names={'FL3', 'NPM1', 'NPM1','CEBPA','KIT','N-RAS',...
    'MLL','MLL', 'MLL','MLL','MLL',...
    'WT1'};
imp_idx=false(1,ngen);
% missing IDH1/2, TET2, DNMT3A, and ASXL1
for i=1:length(imp_genes)
    idx=find(strcmp(genes,imp_genes{i}));
    imp_idx(idx)=1;
end

exp=strsplit(data.textdata{32});
exp=exp(2:end);
nexp=length(exp);
X=data.data;
X=X';

%%
% figure(1);clf;
% hist(X(:), 250);

% gg=sort(X);
% figure(2);clf
% plot(gg);


[ee,idee]=sort(X');
% figure(3);clf
% plot(ee);

[res,perminv]=sort(idee);
q_exp=repmat((1:ngen)'/(ngen+1),1,nexp); %same size as ee
X=q_exp(perminv)';

% keyboard;

% [ee,idee]=sort(X');
% figure(4);clf
% plot(ee);


Xg=norminv(X,0,1);

% gg=sort(X2);
% figure(5);clf
% plot(gg);

% keyboard;

% X2=X;
Xg=X;

%% select genes
% [truefalse, index] = ismember('"GSM332884"', exp) ;

%% remove important genes
% imp_idx=false(1,ngen);
perm=[find(imp_idx) find(~imp_idx) ];
% X2=Xg(:,find(~imp_idx));
X2=Xg(:,perm);
genes=genes(perm);
imp_idx=imp_idx(perm);

%% select genes with more variance
N=500;
N2=N-sum(imp_idx);
S=cov(X2(:,sum(imp_idx)+1:end));
[res ord]=sort(diag(S),'descend');
X3=[X2(:,(sum(imp_idx)+ord(1:N2)))];
X3=[X2(:,1:sum(imp_idx)), X3];

% idx=kmeans(X3', 10);
% [res,order]=sort(idx);
% X4=X3(:,order);

Z = linkage(X3','ward');
% keyboard;
% [H,T,OUTPERM] = dendrogram(Z) ;
[Cres,I]=order_of_tree(Z);
X4=X3(:,I);
genesI=genes(I);
imp_idxI=imp_idx(I);

%% seubseet of 250 genes
X4=X4(:,1:N);

S=cov(X4);

figure(3);clf
subplot(1,2,1);
imagesc(abs(S));
axis square
subplot(1,2,2);
imagesc(imp_idxI');
axis square
colormap jet

X5=[X4(:,imp_idxI) X4(:,~imp_idxI)];
S2=cov(X5);
figure(4);clf
imagesc(abs(S2));
axis square
colormap jet

% figure(4)
% imagesc(Cres(I,:));
% colormap default

S2=cov(X5);
CC=corr(X5);
figure(5);clf
imagesc(abs(CC));
axis square
colormap jet