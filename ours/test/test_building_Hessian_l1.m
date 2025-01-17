
clc;clear all:

p=3;
S=randn(p);
S=S*S';

E=speye(p*p);
pairs=fullfact([p p]);
ii=pairs(:,1)>=pairs(:,2);
pairs=[pairs(ii,1) pairs(ii,2)];

J=1:(2*length(pairs)-p);
I=zeros(1,2*length(pairs)-p);
val=zeros(1,2*length(pairs)-p);

count=1;
for j=1:length(pairs);
    if pairs(j,1)==pairs(j,2)
        J(count)=j;
        I(count)= (pairs(j,2)-1)*p+pairs(j,1);
        val(count)=1;
        count=count+1;
    else
        J(count)=j;
        J(count+1)=j;
        I(count)= (pairs(j,2)-1)*p+pairs(j,1);
        val(count)=1;
        I(count+1)= (pairs(j,1)-1)*p+pairs(j,2);
        val(count+1)=1;
        count=count+2;
    end
end

atoms=sparse(I,J,val);
atoms_l1_sym=[atoms -atoms];
nb_atoms=size(atoms_l1_sym,2);

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

param.max_iter=50;
param.epsilon=1e-3;
param.debug_mode=1;
param.ws=false;

mu=5;
c0=zeros(nb_atoms,1);
new_atom_added=true;
f=q+mu*ones(nb_atoms,1);
[c,A,nbpivot]=asqp(Q,-f,c0,param,new_atom_added);
