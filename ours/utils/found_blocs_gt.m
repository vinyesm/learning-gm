function [ found_blocs ] = found_blocs_gt(ActiveSet,lsupp )
%lsupp ground truth indexes of blocs

nb=size(lsupp,1);
found_blocs=false(1,nb);
for i=1:nb
    currI=lsupp(i,:)';
    found_blocs(i)= isInCell(currI,ActiveSet.I,cell2mat(ActiveSet.k)) ;
end

