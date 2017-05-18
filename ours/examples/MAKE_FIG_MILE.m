
% run('pp_MILE.m');
p=size(Z,1);

%% 
if ~isempty(ActiveSet.alpha)
    Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
else
    Uso=zeros(p,1);
    Dfin=Z1;
end

%% reorder

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
imagesc(abs(-Z1)>1e-6);
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
imagesc(abs(-Z1(I,I))>1e-6);
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
