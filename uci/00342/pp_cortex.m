
n = 77;
A = readtable('Data_Cortex_Nuclear.xls','ReadRowNames',true);

proteins = A.Properties.VariableNames(1:77);
X = A{:,1:77};
set = isnan(X);
Xmean = repmat(nanmean(X,2),1,77);
X(set) = Xmean(set);

C = corr(X);

Z1 = linkage(C,'ward');
I1 = optimalleaforder(Z1,pdist(C));
figure(1);clf;
imagesc(C(I1,I1)); colormap jet; axis square

proteins1 = proteins(I1);
aa =strfind(proteins1,'CaNA_N');
idx = find(not(cellfun('isempty', aa)));

save('cortex', 'C', 'proteins','I1')