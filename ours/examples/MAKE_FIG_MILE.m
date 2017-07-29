
% run('pp_MILE.m');
p=size(Z1,1);

%% 
if ~isempty(ActiveSet.alpha)
    Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
else
    Uso=zeros(p,1);
    Dfin=Z1;
end

%% reorder
% Uso=Uso(:,size(Uso,2):-1:1);
Uso=Uso(:,[1 5 6 7]);
Z2=Uso*Uso';
[I]=grayorder(Uso~=0);

% Z = linkage(full(Uso),'ward');
% [Cres,I]=order_of_tree(Z);
% 
genesI0=genesI;
genesI=genesI(I);
UsoI=Uso(I,:);
Z2II=Z2(I,I);
imp_idxII=imp_idxI(I);
%%

figure(1);clf;
subplot(2,3,1);
imagesc(abs(Uso)>1e-10);
axis square;
subplot(2,3,2);
imagesc(abs(Z2));
title('L')
axis square;
subplot(2,3,3);
imagesc(imp_idxI');
axis square;
subplot(2,3,4);
imagesc(abs(UsoI)>1e-10);
axis square;
subplot(2,3,5);
imagesc(abs(Z2II));
title('L')
axis square;
subplot(2,3,6);
imagesc(imp_idxII');
colormap hot
axis square;

figure(4);clf;
subplot(2,3,1);
imagesc(abs(-Z1)>1e-3);
title('S');
axis square
subplot(2,3,2);
imagesc(abs(Z2));colormap hot
title('L');
axis square
subplot(2,3,3);
imagesc(imp_idxI');
axis square;
subplot(2,3,4);
imagesc(abs(-Z1(I,I))>1e-3);
title('S');
axis square
subplot(2,3,5);
imagesc(abs(Z2(I,I)));colormap hot
title('L');
axis square
subplot(2,3,6);
imagesc(imp_idxII');
colormap hot
axis square;


figure(5);clf;
subplot(2,3,1);
imagesc(abs(Ssl)>1e-3);
title('Ssl');
axis square
subplot(2,3,2);
imagesc(abs(Lsl));colormap hot
title('Lsl');
axis square
subplot(2,3,3);
imagesc(imp_idxI');
axis square;
subplot(2,3,4);
imagesc(abs(Ssl(I,I))>1e-3);
title('Ssl');
axis square
subplot(2,3,5);
imagesc(abs(Lsl(I,I)));colormap hot
title('Lsl');
axis square
subplot(2,3,6);
imagesc(imp_idxII');
colormap hot
axis square

figure(10);clf;
imagesc(-(abs(-Z1(I,I)+eye(p))>1e-3));
colormap hot
axis square;
mkdir('fig');
print('fig/MILE_Som','-depsc')

figure(12);clf;
imagesc(-min(abs(Z2II),.1));
colormap hot
axis square;
mkdir('fig');
set(gca, 'CLim', [-.18, 0]);
print('fig/MILE_Lom','-depsc')


figure(13);clf;
imagesc(-(abs(UsoI)>1e-10));
colormap hot
axis square;
% axis off
pbaspect([1 4 1]);
mkdir('fig');
print('fig/MILE_blocks','-depsc')

figure(15);clf;
imagesc(-abs(Lsl));
colormap hot
axis square;
mkdir('fig');
print('fig/MILE_Lsl_not_ordered','-depsc')

figure(14);clf;
imagesc(-abs(Lsl(I,I)));
colormap hot
axis square;
mkdir('fig');
print('fig/MILE_Lsl_ordered','-depsc')

figure(17);clf;
imagesc(-(abs(Ssl(I,I))>1e-3));
colormap hot
axis square;
mkdir('fig');
print('fig/MILE_Ssl_ordered','-depsc')

% figure(18);
% set(gca, 'CLim', [0 .2]);
% colormap hot;
% colorbar


%%
ActiveSet.alpha=ActiveSet.alpha([1 5 6 7]);
ActiveSet.atoms=ActiveSet.atoms(:,[1 5 6 7]);
save('FIG_MILE','Z1','ActiveSet','genesI','Ssl','Lsl');