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
end
indvalues = zeros(1,length(industry));
for i=1:length(industry)
    indvalues(i)=mapObj(industry{i});
end

%%
% Xn = quantilenorm(X);
Xn0 = log10(X);

Sigma0 = corr(Xn0');
Z0 = linkage(Sigma0,'ward');
I0 = optimalleaforder(Z0,pdist(Sigma0));
figure(1);clf;
imagesc(Sigma0(I0,I0)); colormap jet;


%% remove principal component
% [V,D]= eig(Sigma);
% vi=V(I,1);
% Sigma2 = Sigma(I,I)-D(1,1)*(vi*vi');
% 
% Z = linkage(Sigma2,'ward');
% I2 = optimalleaforder(Z,pdist(Sigma2));
% 
% figure(3);clf;
% imagesc(abs(Sigma2(I2,I2))); colormap jet;
%% Per industry ordering
% 
% K=[];
% for j=1:length(keySet)
%     I3 = find(indvalues==j);
%     Z3 = linkage(Sigma2(I3,I3),'ward');
%     I4 = optimalleaforder(Z3,pdist(Sigma2(I3,I3)));
%     K = [K I3(I4)];
% end

%% extracting block

Xn1 = Xn0(I0(320:end),:);

Sigma1 = corr(Xn1');
Z1 = linkage(Sigma1,'ward');
I1 = optimalleaforder(Z1,pdist(Sigma1));
figure(1);clf;
imagesc(Sigma1(I1,I1)); colormap jet; axis square

Xn = Xn1;
Sigma = Sigma1;
I = I1;

%% splitting train/test

[p,n]=size(Xn);
% nt0 = ceil(2/4*n);
% nt1 = ceil(3/4*n);
% nt2 = ceil(4/4*n);

nb = 20;
J = 1:n;
J = ceil(J/n*nb);
J = (mod(J,2)==1);


Xtrain = Xn(:,J);
Xtest = Xn(:,~J);
Sigma_train = corr(Xtrain');
Sigma_test = corr(Xtest');
figure(50)
subplot(1,2,1);
imagesc(abs(Sigma_train(I,I)));axis square
subplot(1,2,2);
imagesc(abs(Sigma_test(I,I)));axis square

%%
[V,D]=eigs(Sigma_train);
d=diag(D);
[V2,D2]=eigs(Sigma_test);
d2=diag(D2);

figure(30);clf;
plot(abs(V(I,1)*sqrt(d(1)))); hold on
plot(sqrt(V(I,2:3).^2*d(2:3)),'r'); 
hold off

figure(31);clf;
plot(abs(V2(I,1)*sqrt(d2(1)))); hold on
plot(sqrt(V2(I,2:3).^2*d2(2:3)),'r');
hold off

% figure(51)
% subplot(1,2,1);
% imagesc(abs(Sigma_train(I,I)));axis square
% subplot(1,2,2);
% imagesc(abs(Sigma_test(I,I)));axis square

% [p,n]=size(Xn);
% % nt0 = ceil(2/4*n);
% % nt1 = ceil(3/4*n);
% % nt2 = ceil(4/4*n);
% 
% 
% Xtrain = Xn(:,nt0:nt1);
% Xtest = Xn(:,nt1+1:nt2);
% Sigma_train = corr(Xtrain');
% Sigma_test = corr(Xtest');
% figure(50)
% subplot(1,2,1);
% imagesc(abs(Sigma_train(I,I)));axis square
% subplot(1,2,2);
% imagesc(abs(Sigma_test(I,I)));axis square

%%
Sigma=0.5*(Sigma+Sigma');
save('sandp500','Sigma','Sigma_train', 'Sigma_test','names','industry', 'indvalues', 'I');



