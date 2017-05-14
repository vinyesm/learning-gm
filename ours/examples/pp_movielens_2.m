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

addpath ../../ml-latest-small/
addpath ../utils/

%%
ratings=importdata('ratings.csv');
R=ratings.data;
movies=importdata('movies.csv','|');
movies=movies(2:end);
idmovies=cellfun(@getid,movies);

%%
nu=max(R(:,1));
nm=max(idmovies);
idxmovie=1:length(idmovies);

mapping_movies=sparse(idmovies,ones(1,length(idmovies)),idxmovie);

%% all rating mat
Rs=sparse(R(:,1),R(:,2),R(:,3),nu,nm);

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
genre={'Action','Children','Horror'};

rates=sparse(nm,3);
for j=1:3
%     select=(genres(:,idx(j))==1 & sum(Rs(idusers,:)>2)'>0);
    select=~cellfun(@isempty,strfind(movies,genre{j}));
    select=(select & sum(Rs(idusers,idmovies)>2)'>5);
    select= select & ~(idmovies==70);
%     keyboard;
    select=idmovies(select);
%     keyboard;
    for i=select
        rates(i,j)=sum(Rs(idusers,i))/sum(Rs(idusers,i)>0);
    end
end

im=[];
msubset=[];
for i=1:3
    [aa,id]=sort(rates(:,i),'descend');
%     res=idmovies(id);
    im=[im id(1:20)];
    msubset=[msubset mapping_movies(im(:,i))];
end
msubset=full(msubset);
msubset_stack=msubset(:);
im_stack=im(:);

%movies(msubset_stack(43))

%%
X=full(Rs(idusers,im_stack));

% figure(1);clf;
% imagesc(X'); 
% figure(3);clf;
% imagesc(cov(X));
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
inter=real(X3~=0)'*real(X3~=0);
inter=max(inter,1);
S=(X3'*X3);
S=S./inter; %then not PSD
S=.5*(S+S');
S=real(S);
% keyboard;
[U,Ds] = eig(S);
S=U*(Ds.*(Ds>0))*U';
S=.5*(S+S');
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
imagesc(S);
pbaspect([1 1 1]);
title('covariance');
subplot(1,2,2)
imagesc(D);
pbaspect([1 1 1]);
title('inverse covariance');
colormap parula


% %%
% %C=full(cov(X));
% C=cov(X2);
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

