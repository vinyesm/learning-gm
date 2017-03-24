% Toy example
% one latent variable 50%, on grid

clear all; clc;

d=6;
pl=1;
po=d^2;
n=200*po;


pt=pl+po;
density=.5;
pb= round(po*density);
grid=eye(d); %grid of correlations
corr=eye(po);
c=.15;

for i=1:d
    for j=1:d
        idx = (j-1)*d+i;
        if i>1
            corr(idx,idx-1)=c;
            corr(idx-1,idx)=c;
        end
        if i<d
            corr(idx,idx+1)=c;
            corr(idx+1,idx)=c;
        end
        if j>1
            corr(idx,idx-d)=c;
            corr(idx-d,idx)=c;
        end
        if j<d
            corr(idx,idx+d)=c;
            corr(idx+d,idx)=c;
        end
    end
end

% figure(1); clf;
% imagesc(corr);
% pbaspect([1 1 1]);

corr_all= zeros(1+po,1+po);
corr_all(1,1)=1;
corr_all(2:end,2:end)=corr;
corr_all(2:pb+1,1)=ones(pb,1)/sqrt(pb);
corr_all(1,2:pb+1)=ones(pb,1)'/sqrt(pb);
figure(1); clf;
imagesc(corr_all);
pbaspect([1 1 1]);

%%
Dfull=corr_all;
Doo=corr;
Dol=corr_all(2:end,1);
Dmargo=Doo-Dol*Dol';

%%
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


 

