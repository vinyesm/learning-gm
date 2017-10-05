clc; clear all; close all;

addpath DREAM5_NetworkInferenceChallenge_AlternativeDataFormats/net1/
%addpath Network_predictions/Community' 'integration
addpath Network_predictions/Network_predictions/Community' 'integration/
%addpath ../../reorder/
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
% save('adjacency','W');
load('adjacency.mat');

%% cutoff the edges
cut = .713; %1684
A = 1.0*(W > cut);

%% remove isolated nodes
IA = any(A~=0,1);
A = A(IA,IA);

%% Spectral clustering
d = sum(A);
L = diag(d) - A;
Ln = diag(d.^(-0.5))*L*diag(d.^(-0.5));
[V,D]= eig(full(Ln));
ee = diag(D);
idx = find(ee>1e-6,1);

figure(30);
imagesc(max(min(V(:,1:idx), 0.05), -0.05));

%%
% remove big component
 % largest component index
K = 1:idx-1;
ic = kmeans(V(:,K),length(K));

[ics, JJ] = sort(ic);
figure(31);
imagesc(max(min(V(JJ,K), 0.05), -0.05));

figure(32);
plot(max(min(V(JJ,13), 0.05), -0.05),'.');

i0 = 13;
jc = V(:,i0)<-0.01;

% KK = [1:12, 14:idx-1];
% V1 = V(~jc,KK);
% ic1 = kmeans(V1,length(K)-1);
% [ics1, JJ1] = sort(ic1);
% figure(33);
% imagesc(max(min(V1(JJ1,:), 0.05), -0.05));

%%
I = find(~jc);
A1 =  A(~jc,~jc);
p = size(A1,1);
% check 
% norm(A(jc,~jc),'fro')
B = inv(eye(p)-A1./(1.9*norm(full(A1))));
B = B - eye(p);
B = B>1e-6;

figure(34);
imagesc(B);

[C,ia,icc]=unique(B, 'rows');

figure(35);
imagesc(C);

[~,perm] = sort(icc);
Ipc = I(perm);
figure(36);
imagesc(A(Ipc,Ipc)); colormap gray


%% spectral clustering large connex component
I = find(jc);
A2 = A(jc,jc);

d = sum(A2);
L = diag(d) - A2;
Ln = diag(d.^(-0.5))*L*diag(d.^(-0.5));
[V,D]= eig(full(Ln));
ee = diag(D);
idx = find(ee>1e-6,1);
vi = V(:,idx);

figure(37);
plot(sort(vi));

[~,perm] = sort(vi);
nmax = sum(vi < 0);
Igc = I(perm);
Igc = Igc(nmax+1:end);
% Igc = Igc(1:nmax);

%% prendre genes Igc & Ipc
Ic = [Ipc;Igc];
Ac = A(Ic,Ic);
figure(38);
imagesc(Ac); colormap gray

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
expr = BBB.data;
expr = expr(Ic,:);
X  = expr;
%X = quantilenorm(X);
Sigma = cov(X');
Cor = corr(X');

Z = linkage(Cor,'ward');
[Cres,I1]=order_of_tree(Z);
% leafOrder = optimalleaforder(tree,D)
% keyboard;
%[H,T,OUTPERM1] = dendrogram(Z, 0) ;

figure(3);clf;
subplot(1,2,1)
imagesc(abs(Sigma))
axis square
subplot(1,2,2)
imagesc(abs(Sigma(I1,I1)));
axis square

figure(40);clf;
subplot(1,2,1)
imagesc(abs(Cor))
axis square
subplot(1,2,2)
imagesc(abs(Cor(I1,I1)));
axis square

