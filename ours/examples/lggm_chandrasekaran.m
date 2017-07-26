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

%% settings

n0=1000; % number of samples

pl=3; % number of latept variables
po=15;% number of observed variables
pt=po+pl;

slo=.2; % level of sparsity on latent-observed block of conceptration mat
soo=.2; % level of sparsity on latent-observed block of conceptration mat
Too=0*tril(rand(po,po)<soo,-1);
Doo=.2*(Too+Too')+eye(po);
Dol=zeros(po,pl);
Dol(1:5,1)=ones(5,1);
Dol(6:10,2)=ones(5,1);
Dol(11:15,3)=ones(5,1);
Dol=.2*Dol;

%% construction of the concentration matrix
Dfull=zeros(pt);
Dfull(1:pl,1:pl)=eye(pl);
Dfull((pl+1):pt,(pl+1):pt)=Doo;
Dfull(1:pl,(pl+1):pt)=Dol';
Dfull((pl+1):pt,1:pl)=Dol;

Dfull=.5*(Dfull+Dfull');
Dmargo=Doo-Dol*Dol';


descr = {'Plot of the function:';
    'y = A{\ite}^{-\alpha{\itt}}';
    ' ';
    'With the values:';
    'A = 0.25';
    '\alpha = .005';
    't = 0:1000'};

descr = {'Plot of the groundtruth concentration matrix';
    'of the complete model'};

figure(1);clf;
subplot(1,2,1)
axis off;
text(0.5,.5,descr)
pbaspect([1 1 1]);
subplot(1,2,2)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');



if eigs(Dfull,1,'sa')<0
    error('The conceptration matrix of the complete model is not PSD \n')
end

%% sampling data
mu=zeros(1,pt); % vector of means
Xfull=mvnrnd(mu, inv(Dfull), n0)';
X=Xfull((pl+1):pt,:);
S=cov(X');
%S=inv(Dmargo);
S=.5*(S+S');

%% plotting

figure(2);clf;
subplot(2,2,1);
imagesc(abs(Dmargo));
pbaspect([1 1 1]);
title('true marginal conc. mat.');
subplot(2,2,2);
imagesc(abs(inv(S)));
pbaspect([1 1 1]);
title('observed conc. mat.');
subplot(2,2,3);
imagesc(abs(inv(Dmargo)));
pbaspect([1 1 1]);
title('true cov');
subplot(2,2,4);
imagesc(abs(S));
pbaspect([1 1 1]);
title('observed cov');


%% LVGGM Chandrasekaran S-L, (Sparse-Low Rank)

%%
%% set up SDP data in SDPT3 format
%%
Sigma=S;
n=po;
rho=1e-1;
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
beta = 1e-2;
blk{2,1} = 's'; blk{2,2} = n;
At{2,1}  = Atmp;
C{2,1}   = beta*speye(n,n);
%%
blk{3,1} = 'l'; blk{3,2} = 2*n2;
Identity = speye(n2);
At{3,1} = [-Identity,Identity]';
idx = find(Icomp == Jcomp);
ee  = sqrt(2)*ones(m2,1);
if ~isempty(idx); ee(idx) = ones(length(idx),1); end
C{3,1} = rho*[ee; ee];
fprintf('\n Set up data time = %3.2f',etime(clock,ttime));
keyboard;
runPPA = 1;
if (runPPA)
    OPTIONS.smoothing  = 1;
    OPTIONS.scale_data = 0; %% or 2;
    OPTIONS.plotyes    = 0;
    OPTIONS.tol        = 1e-10;
    mu = [1; 0; 0];
    [obj,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,mu,OPTIONS);
    obj = sum(sum(Sigma.*X{1}))-sum(log(eig(X{1})))+rho*sum(sum(abs(X{1})));
    X1 = invD*X{1}*invD; X1 = 0.5*(X1+X1');
    X2 = invD*X{2}*invD; X2 = 0.5*(X2+X2');
end

%% solution
Ssl=X1+X2;
Lsl=X2;

%Usl=chol(Lsl);
[VV DD]=eig(Lsl);
dd=diag(DD);
VV=VV(:,dd>1e-15);
dd=dd(dd>1e-15);
Usl=VV*diag(sqrt(dd));

p=n;
nl=length(dd);

Dsl=zeros(p+nl);
Dsl(1:nl,1:nl)=eye(nl);
Dsl((nl+1):(nl+p),(nl+1):(nl+p))=Ssl;
Dsl(1:nl,(nl+1):(nl+p))=Usl';
Dsl((nl+1):(nl+p),1:nl)=Usl;


figure(3);clf;
subplot(2,2,1);
imagesc(abs(Dmargo));
pbaspect([1 1 1]);
title('true marginal conc. mat.');
subplot(2,2,2);
imagesc(abs(inv(S)));
pbaspect([1 1 1]);
title('observed conc. mat.');
subplot(2,2,3);
imagesc(abs(Ssl+Lsl));
pbaspect([1 1 1]);
title('extimated conc. mat.');
subplot(2,2,4);
axis off
pbaspect([1 1 1]);

descr_gen = {'Formulation S-L where ';
    'with sparse + low rank regularization';
    'solved SD programming';
    ['number of samples n=' num2str(n0)];
    };

descr_par = {['\rho = ' num2str(rho) '  (rho)'];
    ['\beta = ' num2str(beta) '  (beta)'];};

%%
figure(4);clf;
subplot(3,2,1)
axis off;
text(0,.5,descr_gen)
pbaspect([1 1 1]);
subplot(3,2,2)
axis off;
text(0,.5,descr_par)
pbaspect([1 1 1]);
subplot(3,2,3)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
subplot(3,2,4)
imagesc(abs(Dsl));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
subplot(3,2,5)
imagesc(abs(Dfull)>1e-10);
pbaspect([1 1 1]);
title('true support');
subplot(3,2,6)
imagesc(abs(Dsl)>1e-10);
pbaspect([1 1 1]);
title('estimated support');

%filename = ['lggm_chandrasekaran' datestr(datetime('now'),'yyyymmddTHHMMSS') '.ps'];
filename = ['lggm_chandrasekaran' datestr(clock) '.ps'];
print ( '-dpsc2', filename, '-f1' )
print ( '-dpsc2', filename, '-append', '-f2' )
print ( '-dpsc2', filename, '-append', '-f3' )
print ( '-dpsc2', filename, '-append', '-f4' )



