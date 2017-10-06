close all;
addpath individual_stocks_5yr

listing = dir('individual_stocks_5yr');

X = [];
names = {};

for id=1:500%length(listing)
    %A = csvread(listing(3).name); 
    T = readtable(listing(id+3).name);
%     keyboard
    if length(T.Close)==1258
         X = [X;T.Close'];
         aaa = listing(id+3).name;
         nn = strsplit(aaa,'_');
         names = [names nn(1)];
    end
%     X(id,:)=T.Close;
end

%%
%Xn = quantilenorm(X);
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
% [H,T,OUTPERM] = dendrogram(Z) ;
%[Cres,I]=order_of_tree(Z);
I = optimalleaforder(Z,pdist(Sigma));

figure(2);clf;
imagesc(Sigma(I,I)); colormap jet;

save('sandp500','Sigma','names', 'I');

%%
[V,D]= eig(Sigma);
vi=V(I,1);
Sigma2 = Sigma(I,I)-D(1,1)*(vi*vi');

Z = linkage(Sigma2,'ward');
I2 = optimalleaforder(Z,pdist(Sigma2));

figure(3);clf;
imagesc(abs(Sigma2(I2,I2))); colormap jet;

Sigma=0.5*(Sigma2+Sigma2');
I=I2;
save('sandp500','Sigma','names', 'I');

