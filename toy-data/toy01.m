% Toy example
% with three blocks of size resp 20,10,5

clear all; clc; close all;
%% settings

n=1000; % number of samples

ks=[10 10];


pl=length(ks); % number of latent variables
po=sum(ks);% number of observed variables

pt=po+pl;
sigma2l=1;
sigma2lo=1;
scale=ones(pl,1);
% scale=1./sqrt(ks);

%%
% ks=[3 3];
% po=6;
% pl=2;
Uall=zeros(po+pl);
e=[1 .5];
Uo=conv2(eye(po),e,'same');
%Uo(1,[1 po])=e;
Uall(1:pl,1:pl)=eye(pl);
Uall(pl+1:end,:)=[zeros(size(Uo,1),pl) Uo];
for i=1:pl
    if i>1
        i1=sum(ks(1:(i-1)))+1;
    else
        i1=1;
    end
    i2=i1+ks(i)-1;
    Uall((pl+i1):(pl+i2),i)=-.5*ones(ks(i),1);
end

Dfull=zeros(po+pl);
for i=1:size(Uall,1)
    Dfull=Dfull+Uall(i,:)'*Uall(i,:);
end
%%
Doo=Dfull(pl+1:pt,pl+1:pt);
Dol=Dfull(pl+1:pt,1:pl);

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


 