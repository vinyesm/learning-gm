% Movielens data.
% Using the Movielens movie rating data set, we choose the rating  scores 
% given  by  the  most  active  600  users  and  for  the  highest  rated  
% 20 movies from each of the following three genres:
% Horror,Children???s, andAction.
% 
% This results in a 600??60rating matrix with 56percent completeness. 
% We consider the joint distribution of 60 movie rating variables as 
% a LVGGM with three latent variables.  Each user???s rating vector is 
% treated as an i.i.d. sample from the LVGGM. Since the true covariance 
% matrix is unknown, we use the sample covariance matrix
% as a proxy (as n>>p). Each covariance element is weighted by the actual number
% of observations to compensate for the missingness in the data.

clc; clear all

addpath ../../ml-100k/
addpath ../utils/

%%
%19 genres
%1682 movies
%943 users
ng=19;
nm=1682;
nu=943;

%%
info=importdata('u.info');
R=importdata('u.data');
item=importdata('u.item','|');
genres=item.data;
GE=importdata('u.genre');

%% all rating mat
Rs=sparse(R(:,1),R(:,2),R(:,3));

%% 600 most active users
idusers=1:nu;
activity=zeros(nu,1);
for i=1:nu
    activity(i)=sum(R(:,1)==i);
end
[id,aa]=sort(activity,'descend');
idusers=idusers(aa);
idusers=idusers(1:600);


%% selection of movies 
%% highest rated  20 movies from genres
idmovies=1:nm;
idx=[2,5,12]; %2:action 5:children 12:horror

rates=zeros(nm,3);
for j=1:3
    select=(genres(:,idx(j))==1 & sum(Rs(idusers,:)>0)'>1); %rated more than 2 at least 1 times
%     select=(genres(:,idx(j))==1 & sum(Rs(idusers,:)>2)'>5); %rated more than 2 at least 6 times
%     select=(genres(:,idx(j))==1 & sum(Rs(idusers,:)>2)'>10); %(works best)
    select= select & ~(idmovies'==183);
    select= select & ~(idmovies'==234);
    select=idmovies(select);
    for i=select
%         rates(i,j)=sum(Rs(idusers,i))/sum(Rs(idusers,i)>0); %(works best)
%         rates(i,j)=sum(Rs(idusers,i)); %(works best)
        rates(i,j)=sum(Rs(idusers,i)>0); %most rated
    end
end

im=[];
for i=1:3
    [aa,id]=sort(rates(:,i),'descend');
    res=idmovies(id);
    im=[im res(1:20)];
end
idmovies=im;
msubset=item.textdata(idmovies,1:3);
length(unique(im))

% keyboard;

%%
X=full(Rs(idusers,idmovies));

% figure(1);clf;
% imagesc(X');
% 
% figure(3);clf;
% imagesc(cov(X));
% 
% 
% 
% figure(3);clf;
% plot(sum(X>0,2));


%%
X2=X;
% X2(:,42)=[];

% %% weigths (nb of observations)
% w=full(sum(X2>0,2));
% %w=1./(w.^2);
% w=1./w;
% S=weightedcov(X2, w);
% S=full(S);
% D=inv(S);

%% Build covariance

S=zeros(60);
X3=X2;
mu=sum(X2)./sum(X2>0);
X3=X2;
X3(X2==0)=-inf;
X3=bsxfun(@minus,X3,mu);
X3(X3==-inf)=0;
% inter=real(X3~=0)'*real(X3~=0);
% inter=max(inter,1);
inter=ones(60);
S=(X3'*X3);
S=S./inter; %then not PSD
S=.5*(S+S');
S=real(S);
S0=S;
% keyboard;
% [U,Ds] = eig(S);
% S=U*(Ds.*(Ds>0))*U';
% S=.5*(S+S');
% keyboard;

% w=full(sum(X2>0,2));
% w=1./(w.^2);
% w=1./w;
% S=weightedcov(X2, w);

% w=full(sum(X2>0,2));
% S = X3' * (X3 .* repmat(w, 1, 60));                                                  % Weighted Covariance Matrix
% S = 0.5 * (S + S'); 

% w=1./sum(X2>0,2);
% W=diag(w);
% S=(X3'*W*X3)/sum(w);
% S = 0.5 * (S + S'); 

D=inv(S);

figure(3);clf;
subplot(1,2,1)
imagesc(S-diag(diag(S)));
pbaspect([1 1 1]);
title('covariance');
subplot(1,2,2)
imagesc(D);
pbaspect([1 1 1]);
title('inverse covariance');
colormap parula

%%
[U,Ds] = eig(S);
ds=diag(Ds);
U3=U(:,1:3);
d3=ds(1:3);
L3=U3*diag(d3)*U3';
L3=.5*(L3+L3');
figure(4);clf;
subplot(1,1,1)
imagesc(L3);
pbaspect([1 1 1]);
colormap parula

% %%
%C=full(cov(X));
% C=cov(X2);
% % % S=C;
% 
% figure(4);clf;
% subplot(1,2,1)
% imagesc(C);
% pbaspect([1 1 1]);
% title('covariance');
% subplot(1,2,2)
% imagesc(inv(C));
% pbaspect([1 1 1]);
% title('inverse covariance');
% colormap parula
% 
if eigs(S,1,'sa')<0
    fprintf('Covariance not PSD\n');
end

