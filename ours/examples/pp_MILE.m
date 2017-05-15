%% from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13204

addpath ../../MILE/


data=importdata('GSE13204-GPL7473_series_matrix.txt');
td=data.textdata;
% data=importdata('GSE13204_family.soft');
% data=tdfread('GSE13204-GPL7473_series_matrix.txt');

genes=data.textdata(67:1522,1);
exp=strsplit(data.textdata{32});
exp=exp(2:end);
nexp=length(exp);
X=data.data;
X=X';

%%
figure(1);clf;
hist(X(:), 250);

g1=X(1,:);
g1=sort(X);
figure(2);clf
plot(g1);

%%
% X2=log(X);
% 
% figure(2);clf;
% hist(X2(:),250);

%% select genes
[truefalse, index] = ismember('"GSM332884"', exp) ;

%%
S=cov(X);
[res ord]=sort(diag(S),'descend');
X3=X(:,ord(1:500));

% idx=kmeans(X3', 50);
% [res,order]=sort(idx);
% X4=X3(:,order);

Z = linkage(X3','ward');
% [H,T,OUTPERM] = dendrogram(Z) ;
[C,I]=order_of_tree(Z);
X4=X3(:,I);

S=cov(X4);

figure(3);
% imagesc(S(OUTPERM,OUTPERM));
imagesc(S);

figure(4)
imagesc(C(I,:));