function [Hall,fall] = build_Hessian_prox(inputData,param,al1,aom)

debug=0;

Y=inputData.Y;
p=size(Y,1);
mu=param.mu;
lambda=param.lambda;

nb_atoms_l1=size(al1,2);
nb_atoms_om=size(aom,2);
nb_atoms=nb_atoms_l1+nb_atoms_om;

%% building the Hessian
Hall=zeros(nb_atoms,nb_atoms);
fall=zeros(nb_atoms,1);

for i=1:nb_atoms_l1
    Ei=reshape(al1(:,i),p,p);
    if sum(al1(:,i)==2),
        fall(i)=-trace(Y*Ei)+mu*2;
    else
        fall(i)=-trace(Y*Ei)+mu;
    end
    for j=1:i
        Ej=reshape(al1(:,j),p,p);
        Hall(i,j)=trace(Ei*Ej);
        Hall(j,i)=Hall(i,j);
    end
end

for i=1:nb_atoms_om
    Ui=aom(:,i)*aom(:,i)';
    suppi=sum(abs(aom(:,i))>0);
    cf=min(param.cardfun(suppi:end));
%     fall(nb_atoms_l1+i)=-trace(S*Ui)+lambda; %(*)
    fall(nb_atoms_l1+i)=-trace(Y*Ui)+lambda*cf;
    for j=1:i
        Uj=aom(:,j)*aom(:,j)';
        Hall(nb_atoms_l1+i,nb_atoms_l1+j)=trace(Ui*Uj);
        Hall(nb_atoms_l1+j,nb_atoms_l1+i)=Hall(nb_atoms_l1+i,nb_atoms_l1+j);
    end
end

for i=1:nb_atoms_om
    Ui=aom(:,i)*aom(:,i)';
    for j=1:nb_atoms_l1
        Ej=reshape(al1(:,j),p,p);
        Hall(nb_atoms_l1+i,j)=trace(Ui*Ej);
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



