function [Hall,fall] = build_Hessian_l1_SM(inputData,param,al1,aom)
%score matching
debug=0;


S=inputData.X;
p=size(S,1);
mu=param.mu;
lambda=param.lambda;

nb_atoms_l1=size(al1,2);
nb_atoms_om=size(aom,2);
nb_atoms=nb_atoms_l1+nb_atoms_om;

%% building the Hessian
Hall=sparse(nb_atoms,nb_atoms);
fall=sparse(nb_atoms,1);

for i=1:nb_atoms_l1
    Ei=reshape(al1(:,i),p,p);
    if sum(al1(:,i)==2),
%         fall(i)=-trace(S*Ei)+mu*2; %(*) because loss .5*|S^.5(Z1+Z2)S^.5+I|^2
        fall(i)=trace(Ei)+mu*2;
    else
%         fall(i)=-trace(S*Ei)+mu; %(*)
        fall(i)=trace(Ei)+mu;
    end
    for j=1:i
        Ej=reshape(al1(:,j),p,p);
        Hall(i,j)=trace(Ei*S*Ej);
        Hall(j,i)=Hall(i,j);
    end
end

for i=1:nb_atoms_om
    Ui=aom(:,i)*aom(:,i)';
    suppi=sum(abs(aom(:,i))>0);
    cf=min(param.cardfun(suppi:end));
    fall(nb_atoms_l1+i)=trace(Ui)+cf*lambda;
    for j=1:i
        Uj=aom(:,j)*aom(:,j)';
        Hall(nb_atoms_l1+i,nb_atoms_l1+j)=trace(Ui*S*Uj);
        Hall(nb_atoms_l1+j,nb_atoms_l1+i)=Hall(nb_atoms_l1+i,nb_atoms_l1+j);
    end
end

for i=1:nb_atoms_om
    Ui=aom(:,i)*aom(:,i)';
    for j=1:nb_atoms_l1
        Ej=reshape(al1(:,j),p,p);
        Hall(nb_atoms_l1+i,j)=trace(Ui*S*Ej);
        Hall(j,nb_atoms_l1+i)=Hall(nb_atoms_l1+i,j);
    end
end

if debug
    figure(11);clf;
    imagesc(abs(Hall));
    pbaspect([1 1 1]);
    colorbar;
end

end



