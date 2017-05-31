function [Hall,fall] = update_Hessian_l1_sym(S,param,Hold, fold,al1,aom,atom_added)

%atom_added = -1 if none, 1 if l1, 2 if omega

debug=0;


p=size(S,1);
mu=param.mu;
lambda=param.lambda;

nb_atoms_l1=size(al1,2);
nb_atoms_om=size(aom,2);
nb_atoms=nb_atoms_l1+nb_atoms_om;


%% building the Hessian
Hall=sparse(nb_atoms,nb_atoms);
fall=sparse(nb_atoms,1);

if nb_atoms>1
    Hall(1:(nb_atoms-1),1:(nb_atoms-1))=Hold;
    fall(1:(nb_atoms-1))=fold;
end

if atom_added==1
    order=[(1:(nb_atoms_l1-1)) nb_atoms ((nb_atoms_l1):(nb_atoms_l1-1+nb_atoms_om))];
    for i=nb_atoms_l1
        Ei=reshape(al1(:,i),p,p);
        if sum(al1(:,i)==2),
            %         fall(i)=-trace(S*Ei)+mu*2; %(*) because loss .5*|S^.5(Z1+Z2)S^.5+I|^2
            fall(nb_atoms)=+trace(S*Ei)+mu*2;
        else
            %         fall(i)=-trace(S*Ei)+mu; %(*)
            fall(nb_atoms)=+trace(S*Ei)+mu;
        end
        Hall(nb_atoms,nb_atoms)=trace(S*(Ei*(S*Ei)));
        for j=1:(i-1)
            Ej=reshape(al1(:,j),p,p);
            Hall(nb_atoms,j)=trace(S*Ei*S*Ej);
            Hall(j,nb_atoms)=Hall(nb_atoms,j);
        end
        for j=1:nb_atoms_om
            Uj=aom(:,j)*aom(:,j)';
            Hall(nb_atoms_l1-1+j,nb_atoms)=trace(S*(Uj*(S*Ei)));
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
        %     fall(nb_atoms_l1+i)=-trace(S*Ui)+lambda; %(*)
        if param.Sfixed
            fall(nb_atoms)=+trace(S*Ui)+trace(S*param.Sstar*S*Ui)+lambda*cf;
        else
            fall(nb_atoms)=+trace(S*Ui)+lambda*cf;
        end
        for j=1:i
            Uj=aom(:,j)*aom(:,j)';
            Hall(nb_atoms,nb_atoms_l1+j)=trace(S*Ui*S*Uj);
            Hall(nb_atoms_l1+j,nb_atoms)=Hall(nb_atoms,nb_atoms_l1+j);
        end
        for j=1:nb_atoms_l1
            Ej=reshape(al1(:,j),p,p);
            Hall(nb_atoms,j)=trace(S*Ui*S*Ej);
            Hall(j,nb_atoms)=Hall(nb_atoms,j);
        end
    end
end

if debug
    figure(11);clf;
    imagesc(abs(Hall));
    pbaspect([1 1 1]);
    colorbar;
end

end