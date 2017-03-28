pp=188;%192

figure(1);clf
subplot(2,2,1);
imagesc(abs(Dfin1{pp}));
pbaspect([1 1 1])
title(['lat+obs l1+om  best obj lam=' num2str(pair(pp).lambda) ' mu=' num2str(pair(pp).mu)] )
subplot(2,2,2);
imagesc(abs(Dfin2{pp}));
pbaspect([1 1 1])
title(['lat+obs l1+tr best obj  lam=' num2str(pair(pp).lambda) ' mu=' num2str(pair(pp).mu)] )
subplot(2,2,3);
imagesc(abs(Dfin1{pp})>1e-12);
pbaspect([1 1 1])
title(['lat+obs l1+om  best obj lam=' num2str(pair(pp).lambda) ' mu=' num2str(pair(pp).mu)] )
subplot(2,2,4);
imagesc(abs(Dfin2{pp})>1e-12);
pbaspect([1 1 1])
title(['lat+obs l1+tr best obj  lam=' num2str(pair(pp).lambda) ' mu=' num2str(pair(pp).mu)] 