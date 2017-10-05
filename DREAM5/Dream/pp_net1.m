clc; clear all; close all;

addpath DREAM5_NetworkInferenceChallenge_AlternativeDataFormats/net1/
addpath Network_predictions/Community' 'integration
addpath ../../reorder/

CCC=importdata('DREAM5_NetworkInference_Community_Network1.txt');

BBB=importdata('net1_expression_data_avg_t.tsv','\t');

%AAA=xlsread('DREAM5_NetworkInference_Ecoli_Saureus_Modules.xls');

%%
% TF =  unique(CCC.textdata(:,1));
% nedges= length(CCC.data);
% W = sparse([],[],[],1643, 1643,nedges);
% for i=1:nedges
%     str1 = strtok(CCC.textdata(i,1),'G');
%     str2 = strtok(CCC.textdata(i,2),'G');
%     w = CCC.data(i);
%     g1  = str2num(str1{1});
%     g2  = str2num(str2{1});
%     W(g1,g2)=w;
%     W(g2,g1)=w;
% end
% % save('adjacency','W');
load('adjacency.mat');
%%
cut = .713; %1684
A = 1.0*(W > cut);


% %%
% Z = linkage(full(W),'ward');
% % leafOrder = optimalleaforder(tree,D)
% % keyboard;
% [H,T,OUTPERM] = dendrogram(Z, 0) ;
% %[Cres,I]=order_of_tree(Z);
% 
% %%
% figure(1)
% subplot(1,2,1)
% imagesc(A)
% axis square
% subplot(1,2,2)
% imagesc(A(OUTPERM,OUTPERM))
% axis square


%% Spectral clustering
d = sum(A);
L = diag(d) - A;
Ln = diag(d.^(-0.5))*L*diag(d.^(-0.5));
[V,D]= eig(full(Ln));
ee = diag(D);
[val, idx] = min(ee(ee>0));
vi = V(:,idx);
[vsi, perm] = sort(vi);

%%
% ic = kmeans(V(:,idx),2);
% I1 = (ic==1);
% I2 = (ic==2);
% I1 = vi <= median(vi);
% I2 = vi > median(vi);
I1 = vi<=0;
I2 = vi>0;
I=perm;

figure(20);clf;
plot(sort(vi(I1)),'r.');hold on;
plot(sort(vi(I2)),'b.');
%%
figure(21);clf;
imagesc(abs(Ln(I,I))>0);

%%
% G = graph(W(I,I));
% figure(10)
% %plot(G)
% plot(G)
% 
% G1 = graph(W(I1,I1));
% G2 = graph(W(I2,I2));
% figure(11)
% plot(G1)
% figure(12)
% plot(G2)



%%
figure(2)
subplot(1,2,1)
imagesc(W)
axis square
subplot(1,2,2)
imagesc(W(I,I))
axis square

%%
expr = BBB.data;
expr = expr(I1,:);
X = quantilenorm(expr);
Sigma = cov(X');
Cor = corr(X');

Z = linkage(Cor,'ward');
% leafOrder = optimalleaforder(tree,D)
% keyboard;
[H,T,OUTPERM] = dendrogram(Z, 0) ;
%[Cres,I]=order_of_tree(Z);

figure(3);clf;
subplot(1,2,1)
imagesc(abs(Cor))
axis square
subplot(1,2,2)
imagesc(abs(Cor(OUTPERM,OUTPERM)));
axis square
% subplot(2,2,3)
% imagesc(abs(Cor))
% axis square
% subplot(2,2,4)
% imagesc(abs(Cor(I,I)));
% axis square
