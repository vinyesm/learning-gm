addpath ../../ml-100k/

R=importdata('u.data');
IT=importdata('u.item');
GE=importdata('u.genre');

keyboard;
S=cov(F.data);

figure(1);clf;
imagesc(max(abs(S(:)))-abs(S));
colormap bone;