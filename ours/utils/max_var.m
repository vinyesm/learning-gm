function [maxvar]=max_var(Z,ActiveSet,param,inputData )

%
H = gradient(Z,inputData,param);
maxvar=-inf;

for i=1:length(ActiveSet.I)
    %eigenvector associated to largest real eigeinvalue
    currI=ActiveSet.I{i};
    varIJ = norm(H(currI,currI));
    if varIJ>maxvar
        maxvar=varIJ;
    end
end

end
