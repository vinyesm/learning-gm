%% from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13204

addpath ../../MILE/


data=importdata('GSE13204-GPL7473_series_matrix.txt');
td=data.textdata;
% data=importdata('GSE13204_family.soft');
% data=tdfread('GSE13204-GPL7473_series_matrix.txt');

genes=data.textdata(66:1522,1);
ngen=length(genes);
exp=strsplit(data.textdata{32});
exp=exp(2:end);
nexp=length(exp);
X=data.data;
X=X';

%%
% figure(1);clf;
% hist(X(:), 250);

gg=sort(X);
figure(2);clf
plot(gg);


[ee,idee]=sort(X');
figure(3);clf
plot(ee);

[res,perminv]=sort(idee);
q_exp=repmat((1:ngen)'/(ngen+1),1,nexp); %same size as ee
X=q_exp(perminv)';

% keyboard;

[ee,idee]=sort(X');
figure(4);clf
plot(ee);


X2=norminv(X,0,1);

gg=sort(X2);
figure(5);clf
plot(gg);

% keyboard;

X2=X;

%% select genes
[truefalse, index] = ismember('"GSM332884"', exp) ;

%%
S=cov(X2);
[res ord]=sort(diag(S),'descend');
X3=X(:,ord(1:500));

% idx=kmeans(X3', 10);
% [res,order]=sort(idx);
% X4=X3(:,order);

Z = linkage(X3','ward');
% keyboard;
% [H,T,OUTPERM] = dendrogram(Z) ;
[Cres,I]=order_of_tree(Z);
X4=X3(:,I);

%% seubseet of 250 genes
X4=X4(:,1:200);

S=cov(X4);

figure(3);
% imagesc(S(OUTPERM,OUTPERM));
imagesc(abs(S));
colormap jet

% figure(4)
% imagesc(Cres(I,:));
% colormap default