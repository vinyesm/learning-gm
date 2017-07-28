function [ Z,Z1,Z2,Hall,fall, ActiveSet, hist,tau_new] = solve_ps_l1_omega_asqp03( Z,Z1,Z2,ActiveSet,param,inputData,atoms_l1_sym,Hall,fall)
%Using Active Set to solve (PS) problem

debug=1;
debug_update=0;
compute_dg=1;

if param.f==4
Y=inputData.X1*(inputData.X1*Z1*inputData.X1-inputData.Y)*inputData.X1;
end

if param.f==4
    S=sparse(inputData.X1*inputData.X1);
elseif param.f==5
    S=sparse(inputData.X);
end
nbetas=length(ActiveSet.I_l1);

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
new_atom_added=false;
idx_added=-1;
ii=1;
count=1;
cont=true;
H = real(gradient(Z,inputData,param));

while cont
    %% get a new atom ui  [si=ui'*ui] in ActiveSet.atoms  from some C_I I
    % in ActiveSet.I to add to the collection
    % current point is Z
    
    if ActiveSet.atom_count>param.max_nb_atoms
        ActiveSet.atom_count=ActiveSet.atom_count-1;
        cont=false;
        display('maximum number of atoms added');
        %         keyboard;
        break;
    end
    
    if ~isempty(Hall) %~first_pass %|| ActiveSet.atom_count==0
        
        %% Warm start active set
        if new_atom_added
            if new_atom_om
                alpha0=[ActiveSet.beta;ActiveSet.alpha;0];
                idx_added=length(alpha0);
            elseif new_atom_l1
                alpha0=[ActiveSet.beta;0;ActiveSet.alpha];
                idx_added=length(ActiveSet.beta)+1;
            end
        else
            alpha0=[ActiveSet.beta;ActiveSet.alpha];
            idx_added=1;
        end
        
        if ActiveSet.atom_count>param.max_nb_atoms
            fprintf('    nb l1 atoms=%d    nb om atoms=%d\n',length(ActiveSet.I_l1),ActiveSet.atom_count);
            fprintf('    max nb om atoms %d reached \n',param.max_nb_atoms);
            break;
        end
        
        %% active set
        [alph,Jset,npiv]=asqp2(Hall+0*1e-12*speye(length(fall)),-fall,alpha0,param_as,new_atom_added,idx_added);
        
        if debug
            obj0=.5*alpha0'*Hall*alpha0+fall'*alpha0;
            obj1=.5*alph'*Hall*alph+fall'*alph;
            if obj1>obj0
                fprintf('objective increasing in asqp\n');
            end
        end
        
        %% for robustness (too small alpha are deleted)
        
        eps_alph=1e-10;
        
        if(min(abs(alph(Jset))))<eps_alph
            fprintf('small alph\n');
            alph(abs(alph)<eps_alph | alph<0)=0;
            Jset(abs(alph)<eps_alph | alph<0)=0;
        end

        ActiveSet.alpha=alph((nbetas+1):end);
        new_atom_count=sum(Jset);
        ActiveSet.atom_count=new_atom_count;
        ActiveSet.atoms=ActiveSet.atoms(:,Jset);%not necessary (for debbuggging here)
        ActiveSet.alpha=ActiveSet.alpha(Jset);

        Hall=Hall(Jset,Jset);
        fall=fall(Jset);
        
        %% update Z
        Z2=zeros(p);
        nz=find(ActiveSet.alpha>1e-15);
        for j=nz'
            u=ActiveSet.atoms(:,j);
            Z2=Z2+ActiveSet.alpha(j)*(u*u');
        end
        Z=Z1+Z2;
        
        %% Compute objective, loss, penalty and duality gap
        if (param.sloppy==0 || (param.sloppy~=0 && mod(count-1,5)==1)) %&& ~isempty(ActiveSet.alpha)
            if compute_dg
                [loss(ii),pen(ii),obj(ii),dg(ii),time(ii)]=get_val_omega_asqp(Z,Z1,Z2,ActiveSet,inputData,param);
                nb_pivot(ii)=npiv;
                active_var(ii)= sum(ActiveSet.alpha>0);
                dualgap=dg(ii);
                
                if debug && ii>1
                    fprintf(' PS info obj=%f  loss=%f  pen=%f dg=%f  \n', obj(ii),loss(ii), pen(ii), dg(ii));
                    if obj(ii)>obj(ii-1)+param.epsStop
                        fprintf('objective increasing\n');
%                         keyboard;
                    end
                end
            else
                dualgap=inf;
            end
            
            %% Verify sttopping criterion
            H = real(gradient(Z,inputData,param));
            %             [maxIJ,new_row, new_col] = dual_l1_spca(H);
            if param.no_l1
                maxval_l1=-inf;
            else
                [maxval_l1,new_row, new_col] = dual_l1_spca(H);
                cf=1;
            end
            if ~isempty(ActiveSet.I)
                %                 [new_i,new_val,maxvar]=get_new_atom_spca(H,ActiveSet,param);
                [new_i, new_val, maxval_om0]=get_new_atom_spca(H,ActiveSet,param);
                cf=min(param.cardfun(length(new_i):end));
                maxvar=maxval_om0*cf;
            else
                maxvar=0;
            end
         
            if debug
                fprintf('varmax/cf*lambda=%4.2f<1   dg %f<%f \n',maxvar/(cf*param.lambda),dualgap,param.epsStop);
            end            
%             cont=cont && count< param.niterPS;            
%             cont= maxvar/(cf*param.lambda)>1+param.epsStop  && count< param.niterPS;
%             keyboard
            cont=dg(ii)>param.epsStop && count< param.niterPS;
            
            
            if ~cont
%                 keyboard;
            end

            ii=ii+1;         
        end
    else % isempty(Hall)
        if param.no_l1
            maxval_l1=-inf;
        else
            [maxval_l1,new_row, new_col] = dual_l1_spca(H);
        end
        if ~isempty(ActiveSet.I)
            [new_i, new_val, maxval_om0]=get_new_atom_spca(H,ActiveSet,param);
        end
    end
    
    %% get new atom
    
    if cont
        
        new_atom_l1=false;
        new_atom_om=false;
        H = real(gradient(Z,inputData,param));
        if param.no_l1
            maxval_l1=-inf;
        else
            [maxval_l1,new_row, new_col] = dual_l1_spca(H);
            cf=1;
        end

        %% omega atom
        
        cf=1;
        if ~isempty(ActiveSet.I)
            [new_i, new_val, maxval_om0]=get_new_atom_spca(H,ActiveSet,param);
            maxval_om0=maxval_om0*cf;
            anew=sparse(new_i,ones(length(new_i),1),new_val,p,1);
            cf=min(param.cardfun(length(new_i):end));
            if maxval_om0-cf*param.lambda>0
                maxval_om=maxval_om0;
            else
                maxval_om=-inf;
            end
        else
            maxval_om=-inf;
        end
        
        %% no atoms added, break
        if maxval_om==-inf
            fprintf('\n maxval_om is -inf\n');
            break
        end
        
        %% adding new atom
        %         if (new_row==new_col && maxval_l1/param.mu<=maxval_om/(cf*param.lambda)) || (new_row~=new_col && 2*maxval_l1/param.mu<=maxval_om/(cf*param.lambda))
        %% ading omega atom
        if maxval_om<cf*param.lambda
            fprintf('\n not good atom omega d=%f\n',maxval_om);
            keyboard;
        end

        if ActiveSet.atom_count>0
            aom=[ActiveSet.atoms(:,1:ActiveSet.atom_count) anew];
        else
            aom=anew;
        end
        
        if param.f==1
            [Hall_new,fall_new] = update_Hessian_prox(inputData.Y,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
        elseif param.f==4
%             [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
            nat = ActiveSet.atom_count;
            sanew = S*anew;
            supp=sum(abs(anew)>0);
            cf=min(param.cardfun(supp:end));
            Hall_new = zeros(nat+1,nat+1);
            fall_new=zeros(nat+1,1);
            if nat>0
                Hall_new(1:nat,1:nat)=Hall;
                fall_new(1:nat,1)=fall;
                for i=1:nat
                    temp=(sanew'*ActiveSet.atoms(:,i))^2;
                    Hall_new(nat+1,i)=temp;
                    Hall_new(i,nat+1)=temp;
                end                
            end
            Hall_new(nat+1,nat+1)=(sanew'*anew)^2;
            fall_new(nat+1) = anew'*Y*anew +param.lambda*cf;
        elseif param.f==5
            [Hall_new,fall_new] = update_Hessian_l1_SM(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
        end
        %             end
        
        g=Hall_new*[ActiveSet.beta;ActiveSet.alpha;0]+fall_new;
        if g(end)>0
            fprintf('\n Not a descent direction when adding om atom d=%f\n',full(g(end)));
            keyboard;
            break;
        elseif abs(g(end))<1e-10
            fprintf('\n small gradient, stopping\n');
            keyboard;
        else
            ActiveSet.atom_count = ActiveSet.atom_count +1;
            ActiveSet.max_atom_count_reached=max(ActiveSet.max_atom_count_reached,ActiveSet.atom_count);
            %ActiveSet.atomsSupport=[ActiveSet.atomsSupport new_i];
            ActiveSet.atoms(:,ActiveSet.atom_count)=anew;
            Hall=Hall_new;
            fall=fall_new;
            new_atom_added=true;
            new_atom_om=true;
        end
       
        first_pass=false;
        
    end
    
    count=count+1;
end

%redimensioning arrays
% keyboard;
ii=ii-1;
hist.obj = obj(1:ii);
hist.pen =  pen(1:ii);
hist.loss =  loss(1:ii);
hist.dg =  dg(1:ii);
hist.time = time(1:ii);
hist.obj_sup=obj(1);
hist.dg_sup=dg(1);
hist.time_sup=time(1);
hist.nb_pivot=nb_pivot(1:ii);
hist.active_var=active_var(1:ii);

if ii>0
    tau_new=min(param.epsStop,hist.dg(end));
else
    tau_new=param.epsStop;
end

if count>param.niterPS
    fprintf('maximum number of Ps iteration reached, duality gap=%f\n',hist.dg(end));
end


if param.debug
    figure(30);clf;
    plot(hist.obj);
    title('objective');
    keyboard;
end

end


