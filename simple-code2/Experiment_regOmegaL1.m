
clear all
close all
clc

addpath ../ours/TPower_1.0/misc/
addpath ../ours/TPower_1.0/algorithms/TPower/
addpath ../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../spams-matlab-v2.6/build/

run('../DREAM5/Dream/pp_net1');
S05 = sqrt(Sigma);
S05 = (S05+S05')/2;
p= size(Sigma,1);

%%
param.lambda=.02;
param.mu=.002;

%pause()
%keyboard

inputData.X = S05;
inputData.Y = eye(p);
param.k=50;
%set = supp;
param.epsStop=1e-3;
param.maxIter=1000;
param.maxNbAtoms=1000;
param.verbose=2;


[ output, hist ] = regOmegaL1( inputData, param, inf);


figure(1);clf
semilogy(hist.reldgl1,'r');hold on
semilogy(hist.reldgom,'b');
legend('rel dg l1','rel dg om');
title('relative duality gaps of the subproblems');

figure(2);clf
semilogy(hist.reldg,'k');
title('relative global duality gap');

figure(3);clf
semilogy(hist.objective,'k');
title('objective');

figure(17)
imagesc(output.M);
figure(18)
imagesc(M);

figure(19)
imagesc(-full(output.S));
figure(20)
imagesc(S);



%dif=sign(output.S)+sign(S);
%full(output.atoms_u)
