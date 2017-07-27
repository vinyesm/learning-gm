%
% EXPERIMEpt ON non oriepted graph structure learning with latept variables
%
% We build the complete model and sample from it
% We assume latept variables ndependept

%% add paths
clc; clear all; close all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
% addpath('../other');
addpath('../prox');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');
addpath('../../spca_am-master');

HOME = '/Users/marina/Documents/learning-gm/code-from-Kim-Chuan/LogdetPPA-0'; %if my  mac
HOME = '/home/marina/Marina/learning-gm/code-from-Kim-Chuan/LogdetPPA-0';%if lab pc
addpath(strcat(HOME,'/solver/'))
addpath(strcat(HOME,'/solver/mexfun'))
addpath(strcat(HOME,'/util/'))

% addpath LogdetPPA-0/;
% addpath LogdetPPA-0/solver/;
% addpath LogdetPPA-0/solver/mexfun/;
% addpath LogdetPPA-0/util/;
ttime  = clock;

% n=8;
% i=1;
% a=rand(n);
% a=a*a';
% I=[4,5,6];
% k{i}=length(I);
% phi{1}=sparse(I,(1:k{i})',1,n,k{i});
% 
% aII=a(I,I);
% res1=phi{1}*a(I,I)*phi{1}';
% 
% [Iall,Jall] = find(triu(ones(n)));
% tmp = [Iall,Jall];
% m2 = size(tmp,1);
% Icomp = tmp(:,1); Jcomp = tmp(:,2);
% IJ= (Icomp-1)*n+Jall;
% 
% [Ialli,Jalli] = find(triu(ones(k{i})));
% IJi= (Ialli-1)*k{i}+Jalli;
% A=kron(phi{1},phi{1});
% Atriu=A(IJ,IJi);



%% TOY EXAMPLE

% %% data
run('../../toy-data/toy_overlap.m');k=15;

ActiveSet.I={(1:15)',(11:25)',(21:35)',(31:45)'};
k=15;
m=length(ActiveSet.I);


%%
ttime  = clock;
%%
%% set up SDP data in SDPT3 format
%%
Sigma=S;
n=po;
mu=2*.1;
lambda=2*.5;
invD = speye(n,n); 

%
n2 = n*(n+1)/2;
b = zeros(n2,1);
C{1} = Sigma;
blk{1,1} = 's'; blk{1,2} = n;
[Iall,Jall] = find(triu(ones(n)));
tmp = [Iall,Jall];
m2 = size(tmp,1);
Icomp = tmp(:,1); Jcomp = tmp(:,2);
Itmp  = Icomp + Jcomp.*(Jcomp-1)/2;
Atmp  = spconvert([Itmp,[1:m2]',ones(m2,1); n2,m2,0]);
At{1} = Atmp;
%%

[Iall,Jall] = find(triu(ones(n)));
tmp = [Iall,Jall];
m2 = size(tmp,1);
Icomp = tmp(:,1); Jcomp = tmp(:,2);
IJ= (Icomp-1)*n+Jall;

[Ialli,Jalli] = find(triu(ones(k)));
IJi= (Ialli-1)*k+Jalli;

for i=1:m
    blk{1+i,1} = 's'; blk{1+i,2} = k;
    phi=sparse(ActiveSet.I{i},(1:k)',1,n,k);
    A=kron(phi,phi);
    At{1+i,1}  = A(IJ,IJi)';
    C{1+i,1}   = lambda*speye(k,k);
end
%%
blk{m+2,1} = 'l'; blk{m+2,2} = 2*n2;
Identity = speye(n2);
At{m+2,1} = [-Identity,Identity]';
idx = find(Icomp == Jcomp);
ee  = sqrt(2)*ones(m2,1);
if ~isempty(idx); ee(idx) = ones(length(idx),1); end
C{m+2,1} = mu*[ee; ee];
fprintf('\n Set up data time = %3.2f\n',etime(clock,ttime));
runPPA = 1;
if (runPPA)
    OPTIONS.smoothing  = 1;
    OPTIONS.scale_data = 0; %% or 2;
    OPTIONS.plotyes    = 0;
    OPTIONS.tol        = 1e-10;
    eta = [1; zeros(m+1,1)];
    [obj00,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,eta,OPTIONS);
    L=[];
    U=[];
    M=X{1};
    S=M;
    for i=1:m
        I= ActiveSet.I{i};
        L{i}=zeros(n);
        L{i}(I,I)=X{1+i};
        [V,D] = eig(X{1+i});
        d=diag(D);
        ll=sum(d>1e-3);
        u = bsxfun(@times,sqrt(d(d>1e-3))',V(:,d>1e-3));
        U{i}=zeros(n,ll);
        U{i}(I,:)=u;
        S=S+L{i};       
    end
    S(abs(S)<1e-3)=0;
%     obj = sum(sum(Sigma.*X{1}))-sum(log(eig(X{1})))+rho*sum(sum(abs(X{1})));
%     X1 = invD*X{1}*invD; X1 = 0.5*(X1+X1');
%     X2 = invD*X{2}*invD; X2 = 0.5*(X2+X2');
end

%% plot
p=n;
% 
if ~isempty(ActiveSet.I)
    Uso=[];
    for i=1:m
        Uso=[Uso U{i}];
    end
    nl=size(Uso,2);
    Dfin=zeros(p+nl);
    Dfin(1:nl,1:nl)=eye(nl);
    Dfin((nl+1):(nl+p),(nl+1):(nl+p))=S;
    Dfin(1:nl,(nl+1):(nl+p))=Uso';
    Dfin((nl+1):(nl+p),1:nl)=Uso;
else
    Dfin=S;
end

figure(1);clf;
subplot(2,2,1)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
colorbar
subplot(2,2,2)
imagesc(abs(Dfin));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(2,2,3)
imagesc(abs(Dfull)>1e-15);
pbaspect([1 1 1]);
title('true support');
colorbar
subplot(2,2,4)
imagesc(abs(Dfin)>1e-15);
pbaspect([1 1 1]);
title('estimated support');
colorbar
keyboard;
%%
ActiveSet.I={(1:15)',(11:25)',(21:35)',(31:45)'};
ActiveSet.k={15,15,15,15};
% ActiveSet.I={};
% ActiveSet.k={};
param.mu=2*.1;
param.lambda=2*.5;
param.k=15;
param.verbose=1;
param.max_nb_main_loop=100;
param.cardfun=inf*ones(1,p);
param.cardfun(k)=1;
[M,S,L,U,hist,ActiveSet] = logdetPPA_l1_omega(Sigma,param,ActiveSet);

figure(2);clf;
plot(hist.obj_sup);

% 
if ~isempty(ActiveSet.I)
    Uso=[];
    for i=1:length(ActiveSet.I);
        Uso=[Uso U{i}];
    end
    nl=size(Uso,2);
    Dfin2=zeros(p+nl);
    Dfin2(1:nl,1:nl)=eye(nl);
    Dfin2((nl+1):(nl+p),(nl+1):(nl+p))=S;
    Dfin2(1:nl,(nl+1):(nl+p))=Uso';
    Dfin2((nl+1):(nl+p),1:nl)=Uso;
else
    Dfin2=S;
end

figure(3);clf;
subplot(2,2,1)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
colorbar
subplot(2,2,2)
imagesc(abs(Dfin2));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(2,2,3)
imagesc(abs(Dfull)>1e-15);
pbaspect([1 1 1]);
title('true support');
colorbar
subplot(2,2,4)
imagesc(abs(Dfin2)>1e-15);
pbaspect([1 1 1]);
title('estimated support');
colorbar

