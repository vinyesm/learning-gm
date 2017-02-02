function [Hall,fall,Uall] = add_omega_atoms_hessian_l1_sym(inputData,param,H,f,atoms_l1_sym,atoms_om,new_atom)

S=inputData.X1*inputData.X1;
p=size(S,1);
nb_atoms=size(H,2);
nb_atoms_l1=size(atoms_l1_sym,2);

%% building the Hessian
Hall=zeros(nb_atoms+1,nb_atoms+1);

v=zeros(nb_atoms_l1,1);
for i=1:nb_atoms_l1
    Ei=reshape(atoms_l1_sym(:,i),p,p);
    v(i)=new_atom'*(S*Ei*S)*new_atom;
end

unew=real(inputData.X1)*new_atom;
if isempty(atoms_om)
    vnew=[];
else
    vnew=atoms_om'*unew;
    vnew=vnew.*vnew;
end

Uall=[atoms_om unew];
Hall(1:nb_atoms,1:nb_atoms)=H;
Hall(1:nb_atoms,end)=[v;vnew];
Hall(end,1:nb_atoms)=[v;vnew]';
Hall(end,end)=(unew'*unew)^2;
weight=1;
fall=[f;param.lambda*weight-new_atom'*S*unew];

end



