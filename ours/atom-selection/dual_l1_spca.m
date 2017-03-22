function [maxIJ,rowI,colJ] = dual_l1_spca(H)

maxIJ=max(abs(H(:)));
[rowI, colJ] = find(abs(H) == maxIJ);
rowI=rowI(1);
colJ=colJ(1);
if rowI ~= colJ
%     maxIJ=.5*maxIJ; % if outside diagonal it is equivalent to select two atoms at a times
end
end