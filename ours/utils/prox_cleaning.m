function  [Z2,ActiveSet]=prox_cleaning(ActiveSet,thresh)


p=size(ActiveSet.atoms,1);


J=false(1,ActiveSet.atom_count);
for i=1:length(ActiveSet.I)
    supp=ActiveSet.I{i};
    len=length(supp);
    block=zeros(len);
    K=false(1,ActiveSet.atom_count);
    for at=1:ActiveSet.atom_count
        u=ActiveSet.atoms(:,at);
        %         keyboard;
        if len==sum(u~=0) && all(supp==find(u~=0))
            K(at)=1;
            block=block+ActiveSet.alpha(at)*(u(supp)*u(supp)');
        end
    end
    ActiveSet.block{i}=block;
    [U,D]=eig(block);
    D=diag(D);
    KS=abs(D)>thresh;
    KK=find(K);
    KK=KK(1:sum(KS));
    J(KK)=1;
    ActiveSet.atoms(supp,KK)=U(:,KS);
    ActiveSet.alpha(KK)=D(KS);
end
ActiveSet.atom_count=sum(J);
ActiveSet.atoms=ActiveSet.atoms(:,J);
ActiveSet.alpha=ActiveSet.alpha(J);

Z2=zeros(p);
nz=find(ActiveSet.alpha>1e-15);
for j=nz'
    u=ActiveSet.atoms(:,j);
    Z2=Z2+ActiveSet.alpha(j)*(u*u');
end

end
