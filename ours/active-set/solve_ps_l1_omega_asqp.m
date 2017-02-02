function [ Z,Z1,Z2,U,Hall,fall,cardVal, ActiveSet, hist] = solve_ps_l1_omega_asqp( Z,Z1,Z2,ActiveSet,param,inputData,atoms_l1_sym,U,Hall,fall,cardVal)
%Using Active Set to solve (PS) problem

nbetas=size(atoms_l1_sym,2);

param_as.max_iter=1e3;
param_as.epsilon=1e-14;
param_as.debug_mode=false;
param_as.ws=true;

obj = zeros(1,param.niterPS);
pen =  zeros(1,param.niterPS);
loss =  zeros(1,param.niterPS);
dg =  zeros(1,param.niterPS);
time =  zeros(1,param.niterPS);
nb_pivot=zeros(1,param.niterPS);
active_var=zeros(1,param.niterPS);

p=size(Z,1);

if param.debug
    alphaSparsity=[];
    diffalphas=[];
end

new_atom_added=false;
first_pass=1;
i=1;
count=1;
cont=true;

while cont
    %% get a new atom ui  [si=ui'*ui] in ActiveSet.atoms  from some C_I I
    % in ActiveSet.I to add to the collection
    % current point is Z
    
    if ActiveSet.atom_count>param.max_nb_atoms
        display('maximum number of atoms added');
        break;
    end
    
    if ~first_pass || ActiveSet.atom_count==0
        
        if new_atom_added
            alpha0=[ActiveSet.beta;ActiveSet.alpha;0];
        else
            alpha0=[ActiveSet.beta;ActiveSet.alpha];
        end
        
        % active-set
        rcond(Hall+1e-10*eye(length(fall)))
        [alph,Jset,npiv]=asqp(Hall+1e-10*eye(length(fall)),-fall,alpha0,param_as,new_atom_added);
        
        ActiveSet.beta=alph(1:nbetas);
        if length(alph)>nbetas
            Jsetall=Jset;
            Jsetall(1:nbetas)=1;
            Jset=Jset((nbetas+1):end);
            alph=alph((nbetas+1):end);
            new_atom_count=sum(Jset);
            ActiveSet.alpha=alph(Jset);
            ActiveSet.atom_count=new_atom_count;
            ActiveSet.atoms(:,1:new_atom_count)=ActiveSet.atoms(:,Jset);
            U=U(:,Jset);
            Hall=Hall(Jsetall,Jsetall);
            fall=fall(Jsetall);
            cardVal=cardVal(Jset);
            alph
            %             keyboard;
        end
        
        %% Update ActiveSet and Z
        Z1=zeros(p);
        nz=find(ActiveSet.beta>1e-15);
        for j=nz'
            Z1=Z1+ActiveSet.beta(j)*reshape(atoms_l1_sym(:,j),p,p);
        end
        Z2=zeros(p);
        nz=find(ActiveSet.alpha>1e-15);
        for j=nz'
            u=ActiveSet.atoms(:,j);
            Z2=Z2+ActiveSet.alpha(j)*(u*u');
        end
        Z=Z1+Z2;
        
        %% Compute objective, loss, penalty and duality gap
        if (param.sloppy==0 || (param.sloppy~=0 && mod(count,10)==1)) && ~isempty(ActiveSet.alpha)
            [loss(i),pen(i),obj(i),dg(i),time(i)]=get_val_l1_omega_asqp(Z,ActiveSet,inputData,param,cardVal);
            nb_pivot(i)=npiv;
            active_var(i)= sum(ActiveSet.alpha>0);
            cont = (dg(i)>param.PSdualityEpsilon) && count< param.niterPS;
            i=i+1;
        end
    end
    
    if isempty(ActiveSet.I)
        cont=false;
    end
    
    %% get new atom
    if cont
        %% omega atom
        StartY=inputData.Y;
        inputData.Y=StartY-inputData.X1*Z1*inputData.X2;
        [new_i, new_val, maxval_om]=get_new_atom_spca(Z2,ActiveSet,param,inputData);
        inputData.Y=StartY;
        %% l1 sym atom
        StartY=inputData.Y;
        inputData.Y=StartY-inputData.X1*Z2*inputData.X2;
        H = gradient(Z,inputData,param);
        [maxval_l1]=max(abs(H(:)));
        inputData.Y=StartY;
        [new_row new_col] = find(abs(H) == maxval_l1);
        sa=sign(maxval_l1*H(new_row,new_col));
        i1=(new_row-1)*p+new_col;
        i2=(new_col-1)*p+new_row;
        idx_l1 = atoms_l1_sym(i1, :) == sa & atoms_l1_sym(i2, :) == sa;
        keyboard;
        
        if maxval_l1/param.mu>maxval_om/param.lambda
            %adding l1 atom
            if maxval_l1<param.lambda
                fprintf('\n not good atom l1 \n',maxval_l1);
                keyboard;
            end
        else
            %ading omega atom
            if maxval_om<param.lambda
                fprintf('\n not good atom omega d=%f\n',maxval_om);
                keyboard;
            end
        end
        
        %%
        if maxval<param.lambda
            fprintf('\n in solve_ps_spca_asqp Negative directional derivative d=%f\n',maxval);
            keyboard;
        end
        
        anew=sparse(new_i,ones(length(new_i),1),new_val,p,1);
        
        keyboard;
        [Hall_new,fall_new,U_new] = add_omega_atoms_hessian_l1_sym(inputData,param,Hall,fall,atoms_l1_sym,U,anew);
        
        figure(11);clf;
        subplot(1,2,1);
        imagesc(abs(Hall));
        subplot(1,2,2);
        imagesc(abs(Hall_new));
        keyboard;
        
        g=Hall_new*[ActiveSet.beta;ActiveSet.alpha;0]+fall_new;
        
        if g(end)>0
            fprintf('\n in solve_ps_spca_asqp Negative directional derivative d=%f\n',g(end));
            keyboard;
            break;
        else
            ActiveSet.atom_count = ActiveSet.atom_count +1;
            ActiveSet.max_atom_count_reached=max(ActiveSet.max_atom_count_reached,ActiveSet.atom_count);
            ActiveSet.atomsSupport=[ActiveSet.atomsSupport new_i];
            ActiveSet.atoms(:,ActiveSet.atom_count)=anew;
            Hall=Hall_new;
            fall=fall_new;
            U=U_new;
            weight=1;
            cardVal=[cardVal;weight];
            new_atom_added=true;
        end
        
        
        first_pass=false;
        
    end
    
    count=count+1;
end

%redimensioning arrays
i=i-1;
hist.obj = obj(1:i);
hist.pen =  pen(1:i);
hist.loss =  loss(1:i);
hist.dg =  dg(1:i);
hist.time = time(1:i);
hist.obj_sup=obj(1);
hist.dg_sup=dg(1);
hist.time_sup=time(1);
hist.nb_pivot=nb_pivot(1:i);
hist.active_var=active_var(1:i);

if count>param.niterPS
    fprintf('maximum number of Ps iteration reached, duality gap=%f\n',hist.dg(end));
end


if param.debug
    figure(10);clf;
    subplot(1,2,1);
    plot(hist.obj);
    title('objective');
    subplot(1,2,2);
    semilogy(hist.dg);
    title('duality gap');
    keyboard;
end

end


