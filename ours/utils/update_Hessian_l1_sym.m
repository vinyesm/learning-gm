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
%         Ei=reshape(al1(:,i),p,p);
        indx=find(al1(:,i));
        ci=length(indx);
        indx=indx(1);
        [indi,indj]=ind2sub([p p],indx);
        ei=sparse(indi,1,al1(indx,i),p,1);
        ej=sparse(indj,1,1,p,1);
        %         keyboard;
        if sum(al1(:,i)==2),
%             fall(nb_atoms)=+trace(S*Ei)+mu*2;
            fall(nb_atoms)=+ci*(ei'*S*ej)+mu*2;
        else
%             fall(nb_atoms)=+trace(S*Ei)+mu;
            fall(nb_atoms)=+ci*(ei'*S*ej)+mu;
        end
%         Hall(nb_atoms,nb_atoms)=trace(S*(Ei*(S*Ei)));
        Hall(nb_atoms,nb_atoms)=ci^2/2*((ei'*S*ej)^2+(ei'*S*ei)*(ej'*S*ej));
%         if ci==1
%             Hall(nb_atoms,nb_atoms)=(ei'*S*ej)^2;
%         else
%             Hall(nb_atoms,nb_atoms)=2*((ei'*S*ej)^2+(ei'*S*ei)*(ej'*S*ej));
%         end
        for j=1:(i-1)
%             Ej=reshape(al1(:,j),p,p);
            indx=find(al1(:,j));
            cj=length(indx);
            indx=indx(1);
            [indk,indl]=ind2sub([p p],indx);
            ek=sparse(indk,1,al1(indx,j),p,1);
            el=sparse(indl,1,1,p,1);
%             Hall(nb_atoms,j)=trace(S*Ei*S*Ej); 
            Hall(nb_atoms,j)=(ci*cj)/2*((ei'*S*el)*(ej'*S*ek)+(ei'*S*ek)*(ej'*S*el));
            Hall(j,nb_atoms)=Hall(nb_atoms,j);
        end
        for j=1:nb_atoms_om
%             Uj=aom(:,j)*aom(:,j)';
            uj=aom(:,j);
%             Hall(nb_atoms_l1-1+j,nb_atoms)=trace(S*(Uj*(S*Ei)));
            Hall(nb_atoms_l1-1+j,nb_atoms)=ci*(uj'*S*ei)*(uj'*S*ej);
            Hall(nb_atoms,nb_atoms_l1-1+j)=Hall(nb_atoms_l1-1+j,nb_atoms);
%              keyboard;
        end
    end
    Hall=Hall(order,order);
    fall=fall(order);
end

if atom_added==2,
    for i=nb_atoms_om
%         Ui=aom(:,i)*aom(:,i)';
        ui=aom(:,i);
        suppi=sum(abs(aom(:,i))>0);
        cf=min(param.cardfun(suppi:end));
        %     fall(nb_atoms_l1+i)=-trace(S*Ui)+lambda; %(*)
        if param.Sfixed
%             fall(nb_atoms)=+trace(S*Ui)+trace(S*param.Sstar*S*Ui)+lambda*cf;
            fall(nb_atoms)=+(ui'*S*ui)+trace(S*param.Sstar*S*Ui)+lambda*cf;
        else
%             fall(nb_atoms)=+trace(S*Ui)+lambda*cf;
            fall(nb_atoms)=+(ui'*S*ui)+lambda*cf;
        end
        for j=1:i
%             Uj=aom(:,j)*aom(:,j)';
            uj=aom(:,j);
%             Hall(nb_atoms,nb_atoms_l1+j)=trace(S*Ui*S*Uj);
            Hall(nb_atoms,nb_atoms_l1+j)=(ui'*S*uj)^2;
            Hall(nb_atoms_l1+j,nb_atoms)=Hall(nb_atoms,nb_atoms_l1+j);
        end
        for j=1:nb_atoms_l1
%             Ej=reshape(al1(:,j),p,p);
            indx=find(al1(:,j));
            ci=length(indx);
            indx=indx(1);
            [indi,indj]=ind2sub([p p],indx);
            ei=sparse(indi,1,al1(indx,j),p,1);
            ej=sparse(indj,1,1,p,1);
%             Hall(nb_atoms,j)=trace(S*Ui*S*Ej);
            Hall(nb_atoms,j)=ci*(ui'*S*ei)*(ui'*S*ej);
            Hall(j,nb_atoms)=Hall(nb_atoms,j);
%             keyboard;
%             if norm(trace(S*Ui*S*Ej)-ci*(ui'*S*ei)*(ui'*S*ej),'fro')^2>1e-15
%                 fprintf('wrong update\n');
%                 keyboard;
%             end
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