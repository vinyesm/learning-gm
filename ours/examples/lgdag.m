%% add paths
clc; clear all; close all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
addpath('../other');
addpath('../prox');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

load('cov-block.mat')
S=cov;
p=size(cov,1);
lambda=.5;

%% param

param.f=4;
param.diag=0;
param.PSD=true;
param.max_nb_main_loop=30;
param.powerIter=100;
param.stPtPowerIter=1000;
param.niterPS=10000;%5000
param.epsStop=1e-8;
param.PSdualityEpsilon=1e-3;
param.k=0;
param.PSmu=0; %strong convexity
param.verbose=1;
param.debug=0;
param.sloppy=0;
param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
% param.cardfun=inf*ones(1,p);
% param.cardfun(3)=1;
param.cardfun=(1:(p)).^.8;
%%
%S=C;
inputData.X1=S^.5;
inputData.X2=inputData.X1;
inputData.Y=eye(p);
param.lambda=lambda;
%% as quadprog
param.opt='asqp';
[Z_as, ActiveSet_as, hist_as, param_as, flaga_as, it_as] = cgan_spca(inputData,param);
dg_as=hist_as.dg_sup;
tt_as=hist_as.time_sup;
fprintf('lambda=%f\n',lambda);
fprintf('......tt=%f dg=%f\n',hist_as.time(end),hist_as.dg(end));


%% Duality Gap Figure

xplot{1}=tt_as;
yplot{1}=dg_as;
colors = [1 0 0];
legendStr={'cg'};

figure(3);clf
for i=1
loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
hold on
end
title('duality gap lgl');
legend('show','Location','southwest');
grid on
hold off



xplot{1}=hist_as.time;
yplot{1}=hist_as.dg;
colors = [    1 0 0];
legendStr={'cg'};

figure(4);clf
for i=1
loglog(xplot{i},yplot{i},'-','LineWidth',2,'Color',colors(i,:),'DisplayName',legendStr{i}); 
hold on
end
title('duality gap lgl');
legend('show','Location','southwest');
grid on
hold off

figure(5);clf
loglog(tt_as,dg_as,'-','LineWidth',2,'Color',[0 0 0],'DisplayName','dg sup');hold on;
loglog(hist_as.time,hist_as.dg,'-','LineWidth',2,'Color',[1 0 0],'DisplayName','dg sup');hold on;
legend('show','Location','southwest');
grid on
hold off



%%
figure(6);clf
nplots=length(ActiveSet_as.matrix_atoms)+2;
ncol=ceil(sqrt(nplots));
nrow=ceil(nplots/ncol);
iplot=1;
M=zeros(p);
while iplot<nplots-1
subplot(nrow,ncol,iplot);
imagesc(abs(ActiveSet_as.alpha(iplot)*ActiveSet_as.matrix_atoms{iplot}));
M=M+ActiveSet_as.alpha(iplot)*ActiveSet_as.matrix_atoms{iplot};
iplot=iplot+1;
end
subplot(nrow,ncol,iplot);
imagesc(abs(M));
subplot(nrow,ncol,iplot+1);
imagesc(abs(inv(cov)));