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


addpath ../../ml-100k/

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
%% highest rated   20 movies from genres

idx=[2,5,12]; %2:action 5:children 12:horror
idmovies=1:nm;
%select = (sum(genres(:,idx),2)>0);
%select=idmovies(select);
activity=zeros(nm,3);
for j=1:3
    select=(genres(:,idx(j))==1);
    select=idmovies(select);
    for i=select
        activity(i,j)=sum(R(:,2)==i);
    end
end

im=[];
for i=1:3
    [id,aa]=sort(activity(:,i),'descend');
    res=idmovies(aa);
    im=[im res(1:20)];
end
idmovies=im;

%%
Rs=sparse(R(:,1),R(:,2),R(:,3));
X=Rs(idusers,idmovies);
S=full(cov(X));

figure(1);clf;
imagesc(max(abs(S(:)))-abs(S));
colormap bone;