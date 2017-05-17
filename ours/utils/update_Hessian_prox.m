function [Hall,fall] = update_Hessian_prox(Y,param,Hold, fold,al1,aom,atom_added)

%atom_added = -1 if none, 1 if l1, 2 if omega

debug=0;


p=size(Y,1);
mu=param.mu;
lambda=param.lambda;

nb_atoms_l1=size(al1,2);
nb_atoms_om=size(aom,2);
nb_atoms=nb_atoms_l1+nb_atoms_om;


%% building the Hessian
Hall=zeros(nb_atoms,nb_atoms);
fall=zeros(nb_atoms,1);

if nb_atoms>1
    Hall(1:(nb_atoms-1),1:(nb_atoms-1))=Hold;
    fall(1:(nb_atoms-1))=fold;
end

if atom_added==1
    order=[(1:(nb_atoms_l1-1)) nb_atoms ((nb_atoms_l1):(nb_atoms_l1-1+nb_atoms_om))];
    for i=nb_atoms_l1
        Ei=reshape(al1(:,i),p,p);
        if sum(al1(:,i)==2),
            fall(nb_atoms)=-trace(Y*Ei)+mu*2;
        else
            fall(nb_atoms)=-trace(Y*Ei)+mu;
        end
        Hall(nb_atoms,nb_atoms)=trace(Ei*Ei);
        for j=1:(i-1)
            Ej=reshape(al1(:,j),p,p);
            Hall(nb_atoms,j)=trace(Ei*Ej);
            Hall(j,nb_atoms)=Hall(nb_atoms,j);
        end
        for j=1:nb_atoms_om
            Uj=aom(:,j)*aom(:,j)';
            Hall(nb_atoms_l1-1+j,nb_atoms)=trace(Uj*Ei);
            Hall(nb_atoms,nb_atoms_l1-1+j)=Hall(nb_atoms_l1-1+j,nb_atoms);
        end
    end
    Hall=Hall(order,order);
    fall=fall(order);
end

if atom_added==2,
%     keyboard;
    for i=nb_atoms_om
        Ui=aom(:,i)*aom(:,i)';
        suppi=sum(abs(aom(:,i))>0);
        cf=min(param.cardfun(suppi:end));
        fall(nb_atoms)=-trace(Y*Ui)+lambda*cf;
        for j=1:i
            Uj=aom(:,j)*aom(:,j)';
            Hall(nb_atoms,nb_atoms_l1+j)=trace(Ui*Uj);
            Hall(nb_atoms_l1+j,nb_atoms)=Hall(nb_atoms,nb_atoms_l1+j);
        end
        for j=1:nb_atoms_l1
            Ej=reshape(al1(:,j),p,p);
            Hall(nb_atoms,j)=trace(Ui*Ej);
            Hall(j,nb_atoms)=Hall(nb_atoms,j);
        end
    end
end

% keyboard;
if debug
    figure(11);clf;
    imagesc(abs(Hall));
    pbaspect([1 1 1]);
    colorbar;
end

end