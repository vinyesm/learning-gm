clear all; clc; 
load('TOY_OVERLAP.mat');

THRESH=0.02;
Dfin_tr(35,38)=0;
Dfin_tr(38,35)=0;

figure(1);clf;
subplot(3,2,1)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true complete conc. mat.');
colorbar
subplot(3,2,2)
imagesc(abs(Dfull));
pbaspect([1 1 1]);
title('true support');
colorbar
subplot(3,2,3)
imagesc(abs(Dfin).*(abs(Dfin)>THRESH));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(3,2,4)
imagesc(abs(Dfin)>THRESH);
pbaspect([1 1 1]);
title('estimated support');
colorbar
subplot(3,2,5)
imagesc(abs(Dfin_tr).*(abs(Dfin_tr)>THRESH));
pbaspect([1 1 1]);
title('estimated complete conc. mat.');
colorbar
subplot(3,2,6)
imagesc(abs(Dfin_tr)>THRESH);
pbaspect([1 1 1]);
title('estimated support');
colorbar
colormap jet

figure(2);clf;
imagesc(-min(abs(Dfin).*(abs(Dfin)>THRESH),.1));
pbaspect([1 1 1]);
% title('estimated complete conc. mat.');
colormap pink
mkdir('fig');
print('fig/overlap_om','-depsc')

figure(3);clf;
imagesc(-min(abs(Dfin_tr).*(abs(Dfin_tr)>THRESH),.1));
pbaspect([1 1 1]);
% title('estimated complete conc. mat.');
colormap pink
mkdir('fig');
print('fig/overlap_tr','-depsc')

%%
A=ones(size(Dfull,1));
A(Dfull==0)=0;
G = graph(A,'OmitSelfLoops');
figure(4);clf;
plot(G,'Layout', 'circle')
axis equal;
axis off;
mkdir('fig');
print('fig/over_graph','-depsc')
