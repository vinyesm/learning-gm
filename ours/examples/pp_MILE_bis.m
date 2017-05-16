%% from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE13204
clear all; clc;
addpath ../../MILE/


data=importdata('GSE13204-GPL7473_series_matrix.txt');
td=data.textdata;
% data=importdata('GSE13204_family.soft');
% data=tdfread('GSE13204-GPL7473_series_matrix.txt');


%%
genes=data.textdata(66:1522,1);
ngen=length(genes);
exp=strsplit(data.textdata{32});
exp=exp(2:end);
nexp=length(exp);
X0=data.data;

%%
% cf article for classes http://ascopubs.org/doi/pdfdirect/10.1200/JCO.2009.23.4732
leuk_class=strsplit(data.textdata{41},'\t');
leuk_class=leuk_class(2:end);
classes_names=unique(leuk_class);

classes=[];
order=[];
for i=1:length(classes_names)
    class=find(strcmp(leuk_class, classes_names(i)));
    order=[order class];
    classes=[classes i*ones(1,length(class))];
end

X=X0(:,order);

S=cov(X);

figure(3);clf
% imagesc(S(OUTPERM,OUTPERM));
subplot(1,2,1)
imagesc(S);colormap jet;
axis square
subplot(1,2,2)
imagesc(classes');%colormap colorcube;
axis square



% figure(4)
% imagesc(Cres(I,:));
% colormap default