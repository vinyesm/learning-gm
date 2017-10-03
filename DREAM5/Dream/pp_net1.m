clc; clear all; close all;

addpath DREAM5_NetworkInferenceChallenge_AlternativeDataFormats/net1/
addpath Network_predictions/Network_predictions/Community' 'integration
addpath ../reorder/

CCC=importdata('DREAM5_NetworkInference_Community_Network1.txt');

BBB=importdata('net1_expression_data_avg_t.tsv','\t');

%AAA=xlsread('DREAM5_NetworkInference_Ecoli_Saureus_Modules.xls');

%%
TF =  unique(CCC.textdata(:,1));
nedges= length(CCC.data);
adjacency = sparse([],[],[],1643, 1643,nedges);
for i=1:nedges
    str1 = strtok(CCC.textdata(i,1),'G');
    str2 = strtok(CCC.textdata(i,2),'G');
    g1  = str2num(str1{1});
    g2  = str2num(str2{1});
    adjacency(g1,g2)=1;
end

%%
Z = linkage(full(adjacency'),'ward');
% keyboard;
[H,T,OUTPERM] = dendrogram(Z) ;
%[Cres,I]=order_of_tree(Z);

%%
expr = BBB.data;
X = quantilenorm(expr);
Sigma = cov(X');

% figure(1)
% for i=1:10
%     hist(X(:,i),40);
% end