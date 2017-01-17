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


Z=proj_omega(A,0.01);

figure(1);clc;
subplot(1,2,1);
imagesc(abs(A0));
subplot(1,2,2);
imagesc(abs(Z));