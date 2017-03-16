function [maxvar, kmaxvar]=max_var(Z,ActiveSet,param,inputData )

%
H = gradient(Z,inputData,param);
maxvar=-inf;
kmaxvar=0;

% for i=1:length(ActiveSet.I)
%     %eigenvector associated to largest real eigeinvalue
%     currI=ActiveSet.I{i};
%     varIJ = norm(H(currI,currI));
%     if varIJ>maxvar
%         maxvar=varIJ;
%         kmaxvar=i;
%     end
% end

for i=1:length(ActiveSet.I)
    %eigenvector associated to largest real eigeinvalue
    currI=ActiveSet.I{i};
    %varIJ = norm(H(currI,currI));
    A=.5*(H(currI,currI)+H(currI,currI)');
    varIJ = eigs(-A,1,'la');
    if varIJ>maxvar
        maxvar=varIJ;
        kmaxvar=i;
    end
end
% keyboard;

end
