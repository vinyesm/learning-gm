function [new_i,new_val,val]=get_new_atom_spca(H,ActiveSet,param )
% keyboard;
%
% H = gradient(Z,inputData,param);
maxval=-inf;
val=-inf;
new_i=[];
new_val=[];

for i=1:length(ActiveSet.I)
    %eigenvector associated to largest real eigeinvalue
    S=ActiveSet.I{i};
    cf=min(param.cardfun(length(S):end));
    Hs=H(S,S);
    Hs=0.5*(Hs+Hs');
    emin=norm(Hs, 'fro');
%     [v,d,flag]=eigs(-Hs,1,'la');
    [v,d,flag]=eigs(-Hs+emin*speye(length(S)),1,'lm');
    if flag
        fprintf('in get_new_atom_spca eigs has not converged.\n');
%         d=-v'*Hs*v;
%         keyboard;
        d=-1;  
    end
    v=real(v);
    d=real(d);
    d=d-emin;
    if d/cf>maxval
        maxval=d/cf;
        val=d;
        new_i=S;
        new_val=v/sqrt(cf);
    end
end
% fprintf('\n%f',maxval);
if maxval<=0
    fprintf('in get_new_atom_spca : Largest eigenvalue is negative or zero\n');
%     keyboard
end

if isempty(new_i)
    fprintf('in get_new_atom_spca : new_i is empty\n');
%     keyboard
end


% anew=sparse(maxatom);


if ~isreal(v)
    error('selected new atom is not real\n')
end


end
