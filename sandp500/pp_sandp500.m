close all;clear all; clc;
addpath individual_stocks_5yr


listing = dir('individual_stocks_5yr');

X = [];
industry = [];
names = {};

keySet =   {'Consumer Discretionary', 'Consumer Staples','Utilities' ...
            'Energy','Industrials', 'Materials', ...
            'Real Estate', 'Financials',...
            'Information Technology', 'Telecommunication Services',...
            'Health Care'};
valueSet = 1:length(keySet);
mapObj = containers.Map(keySet,valueSet);

%%
TI = readtable('industries.csv');
%%
for id=1:500%length(listing)
    %A = csvread(listing(3).name); 
    T = readtable(listing(id+3).name);
%     keyboard
    
    if length(T.Close)==1258
        aaa = listing(id+3).name;
        nn = strsplit(aaa,'_');
        nn =nn(1);
        ind=find(ismember(TI.Code,nn));
         if ~isempty(ind)
         X = [X;T.Close'];
         names = [names nn];
         industry = [industry TI.Industry(ind)];
         end
    end
%     X(id,:)=T.Close;
end

%%
% Xn = quantilenorm(X);
%Xn = log10(X);
Xn = log10(X);
Sigma = corr(Xn');
d= diag(Sigma);
% 
% figure(5);
% plot(d);
% Ismall = d<5000;

% Sigma=Sigma(Ismall,Ismall);
% names= names(Ismall);

figure(1);clf;
imagesc(Sigma); colormap jet;

%%

% Z = linkage(Sigma,'ward');
% [Cres,I]=order_of_tree(Z);
Z = linkage(Sigma,'ward');
I = optimalleaforder(Z,pdist(Sigma));

figure(2);clf;
imagesc(Sigma(I,I)); colormap jet;

save('sandp500','Sigma','names', 'I');

%% remove principal component
[V,D]= eig(Sigma);
vi=V(I,1);
Sigma2 = Sigma(I,I)-D(1,1)*(vi*vi');

Z = linkage(Sigma2,'ward');
I2 = optimalleaforder(Z,pdist(Sigma2));

figure(3);clf;
imagesc(abs(Sigma2(I2,I2))); colormap jet;

%% Per industry ordering

indvalues = zeros(1,length(industry));
for i=1:length(industry)
    indvalues(i)=mapObj(industry{i});
end
% [~, K]=sort(indvalues);
K=[];
for j=1:length(keySet)
    I3 = find(indvalues==j);
    Z3 = linkage(Sigma2(I3,I3),'ward');
    I4 = optimalleaforder(Z3,pdist(Sigma2(I3,I3)));
    K = [K I3(I4)];
end

figure(4);clf;
subplot(1,2,1)
imagesc(abs(Sigma(K,K))); colormap jet;
subplot(1,2,2)
imagesc(indvalues(K)'); colormap jet;
pbaspect([50 length(industry) 1])

figure(5);clf;
subplot(1,2,1)
imagesc(abs(Sigma2(K,K))); colormap jet;
subplot(1,2,2)
imagesc(indvalues(K)'); colormap jet;
pbaspect([50 length(industry) 1])

figure(6); clf;
subplot(1,2,1)
imagesc(abs(Sigma2(I2,I2))); colormap jet;
subplot(1,2,2)
imagesc(indvalues(I2)'); colormap jet;

figure(7); clf;
subplot(1,2,1)
imagesc(abs(Sigma(I,I))); colormap jet; 
subplot(1,2,2)
imagesc(indvalues(I)'); colormap jet;


%%


%%
[p,n]=size(Xn);
nt = ceil(2/3*n);
Xtrain = Xn(:,1:nt);
Xtest = Xn(:,nt+1:end);
Sigma_train = corr(Xtrain');
Sigma_test = corr(Xtest');
figure(50)
subplot(1,2,1);
imagesc(abs(Sigma_train(I,I)));
subplot(1,2,2);
imagesc(abs(Sigma_test(I,I)));

%%
Sigma=0.5*(Sigma2+Sigma2');
% I=I2;
save('sandp500','Sigma','Sigma_train', 'Sigma_test','names','industry', 'indvalues', 'I','K');



