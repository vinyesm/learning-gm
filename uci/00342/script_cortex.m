
clear all
close all
clc

addpath ../../ours/TPower_1.0/misc/
addpath ../../ours/TPower_1.0/algorithms/TPower/
addpath ../../ours/TPower_1.0/algorithms/PathSPCA/PathSPCA/
addpath ../../spams-matlab-v2.6/build/

addpath ../../simple-code2

load('cortex');

param.lambda=.4;
param.mu=.1;
param.epsStop=1e-5;

%pause()
%keyboard

inputData.X = C;
inputData.Y = eye(77);
param.k=10;
%param.epsStop=1e-3;
param.maxIter=10000;
param.maxNbAtoms=1000;
param.verbose=2;


[ output, hist ] = regOmegaL1( inputData, param, inf );


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
imagesc(abs(output.M(I1,I1)));

figure(19)
imagesc(abs(output.S(I1,I1)));



%dif=sign(output.S)+sign(S);
%full(output.atoms_u)
