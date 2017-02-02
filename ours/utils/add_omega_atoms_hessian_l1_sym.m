function [Qall,qall] = add_omega_atoms_hessian_l1_sym(S,Q,q,H,f,atoms_l1_sym,new_atom)
p=size(S,1);
nb_atoms_l1=size(Q,2);
nb_atoms_om=size(H,2);
nb_atoms=nb_atoms_l1+nb_atoms_om;

%% building the Hessian

Q=zeros(nb_atoms);
q=zeros(nb_atoms,1);
for i=1:nb_atoms
    Ei=reshape(atoms_l1_sym(:,i),p,p);
    q(i)=-trace(S*Ei);
    for j=1:i
        Ej=reshape(atoms_l1_sym(:,j),p,p);
        Q(i,j)=trace(S*Ei*S*Ej);
        Q(j,i)=Q(i,j);
    end
end


end



