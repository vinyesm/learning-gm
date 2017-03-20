function [maxvar, kmaxvar,kmax]=dual_om_k_spca(H,ActiveSet,param)

%
maxvar=-inf;
kmaxvar=0;

for i=1:length(ActiveSet.I)
    %eigenvector associated to largest real eigeinvalue
    currI=ActiveSet.I{i};
    supp=length(currI);
    %varIJ = norm(H(currI,currI));
    if (param.cardfun(supp) ~= inf)
        cf=param.cardfun(supp);
        HI=.5*(H(currI,currI)+H(currI,currI)');
        varIJ = eigs(-HI,1,'la');
        varIJ = varIJ/cf;
        if varIJ>maxvar
            maxvar=varIJ;
            kmaxvar=i;
            kmax=supp;
        end
    end
%     A=.5*(H(currI,currI)+H(currI,currI)');   
end
% keyboard;

end
