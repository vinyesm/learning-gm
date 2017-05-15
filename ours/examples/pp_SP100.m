% SP100 dataset
% sic code : Standard Industrial Classification 
% description in https://www.sec.gov/info/edgar/siccodes.htm

clc; clear all

addpath ../../latentTree/data
addpath ../utils/

load('cov_sp100.mat');

S=cov_mat;

D=inv(S);

figure(3);clf;
subplot(1,2,1)
imagesc(S);
pbaspect([1 1 1]);
title('covariance');
subplot(1,2,2)
imagesc(D);
pbaspect([1 1 1]);
title('inverse covariance');
colormap parula