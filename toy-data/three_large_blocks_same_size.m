% Toy example
% with three blocks of size 5

clear all; clc; close all;
%% settings

n=500; % number of samples

k=10;  %size of the block
B=5;   %nb of blocks

pl=B; % number of latept variables
po=k*B;% number of observed variables

pt=po+pl;
sigma2l=.2;
sigma2lo=.2;
scale=1/sqrt(k);
scale=1;

% slo=.2; % level of sparsity on latent-observed block of conceptration mat
% soo=.2; % level of sparsity on latent-observed block of conceptration mat
Doo=eye(po)/sigma2lo;
Dol=zeros(po,pl);
for i=1:B
    i1=k*(i-1)+1;
    i2=i1+k-1;
    Dol(i1:i2,i)=ones(k,1);
end
Dol=-Dol*scale/sigma2lo;
Dll=eye(pl)*(1/sigma2l+k*scale^2/sigma2lo);

%% construction of the conceptration matrix
Dfull=zeros(pt);
Dfull(1:pl,1:pl)=Dll;
Dfull((pl+1):pt,(pl+1):pt)=Doo;
Dfull(1:pl,(pl+1):pt)=Dol';
Dfull((pl+1):pt,1:pl)=Dol;

Dfull=.5*(Dfull+Dfull');
Dmargo=Doo-Dol*Dol';


descr = {'Plot of the groundtruth '
    'concentration matrix';
    'of the complete model'};

figure(1);clf;
subplot(1,2,1)
axis off;
text(0,.5,descr)
pbaspect([1 1 1]);
subplot(1,2,2)
imagesc(abs(Dfull));
colorbar;
pbaspect([1 1 1]);
title('true complete conc. mat.');

if eigs(Dfull,1,'sa')<0
    error('The conceptration matrix of the complete model is not PSD \n')
end
if eigs(Dmargo,1,'sa')<0
    error('Dmargo matrix of the complete model is not PSD \n')
end

%% sampling data
mu=zeros(1,pt); % vector of means
Xfull=mvnrnd(mu, inv(Dfull), n)';
X=Xfull((pl+1):pt,:);
S=cov(X');
% S=inv(Dmargo);

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


 