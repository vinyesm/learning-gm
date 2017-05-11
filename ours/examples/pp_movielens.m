% Movielens data.
% Using the Movielens movie rating data set, we choose the rating  scores 
% given  by  the  most  active  600  users  and  for  the  highest  rated  
% 20 movies from each of the following three genres:
% Horror,Children’s, andAction.
% 
% This results in a 600×60rating matrix with 56percent completeness. 
% We consider the joint distribution of 60 movie rating variables as 
% a LVGGM with three latent variables.  Each user’s rating vector is 
% treated as an i.i.d. sample from the LVGGM. Since the true covariance 
% matrix is unknown, we use the sample covariance matrix
% as a proxy (as n>>p). Each covariance element is weighted by the actual number
% of observations to compensate for the missingness in the data.


addpath ../../ml-100k/

R=importdata('u.data');
IT=importdata('u.item');
GE=importdata('u.genre');

keyboard;
S=cov(F.data);

figure(1);clf;
imagesc(max(abs(S(:)))-abs(S));
colormap bone;