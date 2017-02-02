function [new_i,new_val,maxval]=get_new_atom_spca(Z,ActiveSet,param,inputData )

%
p=size(Z,1);
H = gradient(Z,inputData,param);
maxval=-inf;

for i=1:length(ActiveSet.I)
    %eigenvector associated to largest real eigeinvalue
    S=ActiveSet.I{i};
    cf=param.cardfun(length(S));
    Hs=H(S,S);
    Hs=0.5*(Hs+Hs');
    [v,d]=eigs(-Hs,1,'la');
    v=real(v);
    d=real(d);
    if d/cf>maxval
        maxval=d/cf;
        new_i=S;
        new_val=v/sqrt(cf);
%         maxatom=zeros(p,1);
%         maxatom(S)=v;
    end
end
% fprintf('\n%f',maxval);
if maxval<=0
    display('in get_new_atom_spca : Largest eigenvalue is negative or zero\n');
    %keyboard
end


% anew=sparse(maxatom);


if ~isreal(v)
    error('selected new atom is not real\n')
end


end
