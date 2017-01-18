%%
%testing projection onto Omega

%% add paths
clc; clear all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../prox');
addpath('../utils');
addpath('../other');

%%
n=40+20+10+5;
A0=zeros(n);
W=randn(n);
W=(W+W')/2;
A0(1:40,1:40)=rand_sym_mat(40,2);
A=A0+0*W;
param.stPtPowerIter=100;
param.powerIter=200;

B=rand_sym_mat(30,5);
A=B;
eig(B)

keyboard;
Zproj=proj_omega(A,0.01);
Zprox=prox_omega(A,0.01);

figure(1);clc;
subplot(1,3,1);
imagesc(abs(A0));
title('X0')
subplot(1,3,2);
imagesc(abs(Zproj));
title('proj')
subplot(1,3,3);
imagesc(abs(Zprox));
title('prox')

%% test soft_threshold
% A=rand(5)-.5;
% B=soft_threshold(A,.2);