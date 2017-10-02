% % Breast Cancer NKI
% Processing the breast cancer gene expression data (in breastCancer.R). 
% The breast cancer data set is maintained by Dana-Farber Cancer Institute 
% at Harvard University, and is available through their R package: 
% breastCancerNKI version 1.3.1 [70]. 
% We removed genes with > 10% missing values. We imputed the remaining 
% missing values using the R package impute (version 1.36.0) [71,
% 111]. 
% We projected the gene expression levels of each gene to the quantiles 
% of a standard normal.
% There were 24,158 genes remaining in the data set after filtering

clc; clear all;
load('Rosetta.mat', 'exprs');
load('Rosetta.mat', 'names');
load('Rosetta.mat', 'erdata');
X = exprs.data;
% Xn = quantilenorm(X,'DISPLAY',true);
Xn0 = quantilenorm(X);
Xn = quantilenorm(Xn0')';
Xn = Xn(1:500,:);
[n,p] = size(Xn);
% C = Xn*Xn';

%%
%% select genes with more variance
% N=250;
% N2=N-sum(imp_idx);
% S=cov(X2(:,sum(imp_idx)+1:end));
% [res ord]=sort(diag(S),'descend');
% X3=[X2(:,(sum(imp_idx)+ord(1:N2)))];
% X3=[X2(:,1:sum(imp_idx)), X3];


%%
Z = linkage(Xn,'ward');
% [H,T,OUTPERM] = dendrogram(Z) ;
[Cres,I]=order_of_tree(Z);
X4=Xn(I,:);

%S=cov(X');
S=cov(X4');

%%
figure(1);clf
% imagesc(S(OUTPERM,OUTPERM));
imagesc(min(abs(S),0.5)); colormap hot;
axis square
Sigma = S;
save('BC','Sigma');

figure(2);clf
% imagesc(S(OUTPERM,OUTPERM));
imagesc(abs(corr(X4'))); colormap hot;
axis square

% x = sort(X(:,1:20));
% %pr = [0.025 0.25 0.50 0.75 0.975];
% pr = (1:100)/101;
% Q = quantile(x,pr);
% N0 = norminv(pr)';

% figure(1);clf;
% subplot(1,2,1)
% plot(x)
% title('scatter points')
% subplot(1,2,2)
% % plot(N0,N0,'k--'); hold on
% plot(Q,N0)
% title('Q-Q plot')

%%
for i=1:80
    fprintf('%s\n',names{i});
end