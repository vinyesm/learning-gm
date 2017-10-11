function [ns, nl] = nbModelParam(S, U)
    p= size(S,1);
    ns = p + sum(sum(abs(S-diag(diag(S)))>.001))/2;
    nl=0;
    for i=1:length(U)
        for j=1:size(U{i},2)
            nl = nl + sum(abs(U{i}(:,j))>.001);
        end
    end
end