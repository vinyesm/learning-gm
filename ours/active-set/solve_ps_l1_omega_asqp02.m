function [ Z,Z1,Z2,Hall,fall, ActiveSet, hist] = solve_ps_l1_omega_asqp02( Z,Z1,Z2,ActiveSet,param,inputData,atoms_l1_sym,Hall,fall)
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
i=1;
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
            obj0=.5*alpha0'*Hall*alpha0+fall'*alpha0
            obj1=.5*alph'*Hall*alph+fall'*alph
            if obj1>obj0
                fprintf('objective increasing in asqp\n');
            end
        end
        
        %% for robustness (too small alpha are deleted)
        
        eps_alph=0;1e-8;
        
        if(min(abs(alph(Jset))))<eps_alph
            fprintf('small alph\n');
            alph(abs(alph)<eps_alph | alph<0)=0;
            Jset(abs(alph)<eps_alph | alph<0)=0;
        end
        
        if nbetas>0
            Jbeta=Jset(1:nbetas);
            ActiveSet.beta=alph(1:nbetas);
            ActiveSet.beta=ActiveSet.beta(Jbeta);
            ActiveSet.I_l1=ActiveSet.I_l1(Jbeta);
        end
        
        if length(alph)>nbetas
            Jalpha=Jset((nbetas+1):end);
            ActiveSet.alpha=alph((nbetas+1):end);
            new_atom_count=sum(Jalpha);
            ActiveSet.atom_count=new_atom_count;
            ActiveSet.atoms=ActiveSet.atoms(:,Jalpha);%not necessary (for debbuggging here)
            ActiveSet.alpha=ActiveSet.alpha(Jalpha);
        end
        
        nbetas=length(ActiveSet.beta);
        Hall=Hall(Jset,Jset);
        fall=fall(Jset);
        
        %% Update ActiveSet and Z
%         if ~param.no_l1
%             if ~param.Sfixed
%                 Z1=zeros(p);
%                 nz=find(ActiveSet.beta>1e-15);
%                 for j=nz'
%                     Z1=Z1+ActiveSet.beta(j)*reshape(atoms_l1_sym(:,ActiveSet.I_l1(j)),p,p);
%                 end
%             end
%         end
        Z2=zeros(p);
        nz=find(ActiveSet.alpha>1e-15);
        for j=nz'
            u=ActiveSet.atoms(:,j);
            Z2=Z2+ActiveSet.alpha(j)*(u*u');
        end
        Z=Z1+Z2;
        
        %% Compute objective, loss, penalty and duality gap
        if (param.sloppy==0 || (param.sloppy~=0 && mod(count,100)==1)) %&& ~isempty(ActiveSet.alpha)
            if compute_dg
                [loss(i),pen(i),obj(i),dg(i),time(i)]=get_val_l1_omega_asqp(Z,ActiveSet,inputData,param);
                nb_pivot(i)=npiv;
                active_var(i)= sum(ActiveSet.alpha>0);
                dualgap=dg(i);
                %cont = (dg(i)>param.PSdualityEpsilon) && count< param.niterPS;
                
                if debug && i>1
                    fprintf(' PS info obj=%f  loss=%f  pen=%f penl1=%f pen_om=%f dg=%f  \n', obj(i),loss(i), pen(i), param.lambda*sum(ActiveSet.alpha),param.mu*sum(ActiveSet.beta), dg(i));
                    if obj(i)>obj(i-1)+param.epsStop
                        fprintf('objective increasing\n');
                        keyboard;
                    end
                end
            else
                dualgap=inf;
            end
            i=i+1;
            
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
            %% extra condition
            if compute_dg
                omega=pen(i-1);
            else
                cfa=zeros(ActiveSet.atom_count,min(ActiveSet.atom_count,1));
                card=sum(ActiveSet.atoms(:,1:ActiveSet.atom_count)~=0);
                for i=1:ActiveSet.atom_count
                    cfa(i)=min(param.cardfun(card(i):end));
                end
                if ~isempty(ActiveSet.alpha)
                    omega = param.lambda*dot(cfa,ActiveSet.alpha);
                else
                    omega = 0;
                end
            end
            
            dotHZ=trace(-H*Z);
            dualomega=max(maxvar/(cf*param.lambda),maxval_l1/param.mu);
            cond=omega*dualomega - dotHZ;
            epscond=param.epsStop;
           
            if cond<epscond
                cont=false;
            end
            
            if debug
                fprintf('   maxIJ/mu=%4.2f<1     varmax/cf*lambda=%4.2f<1   dg/eps=%4.2f<1  cond=%4.2f<%4.2f\n',maxval_l1/param.mu,maxvar/(cf*param.lambda),dualgap/param.PSdualityEpsilon,cond,epscond);
            end
            
            cont=cont && count< param.niterPS;
            
            cont=count<10;
            
            
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
        
        %% omega atom
        
        cf=1;
        if ~isempty(ActiveSet.I)
            %             [new_i, new_val, maxval_om0]=get_new_atom_spca(H,ActiveSet,param);
            maxval_om0=maxval_om0*cf;
            anew=sparse(new_i,ones(length(new_i),1),new_val,p,1);
            cf=min(param.cardfun(length(new_i):end));
            if maxval_om0-cf*param.lambda>param.epsStop/10
                maxval_om=maxval_om0;
            else
                maxval_om=-inf;
            end
        else
            maxval_om=-inf;
        end
        
        %         %% l1 sym atom
        %
        %         %         [maxval_l1,new_row, new_col] = dual_l1_spca(H);
        %
        %         %         if maxval_l1>param.mu*(1+param.epsStop)
        %         if maxval_l1-param.mu>param.epsStop/10
        %             sa=-sign(H(new_row,new_col));
        %             i1=(new_row-1)*p+new_col;
        %             i2=(new_col-1)*p+new_row;
        %             idx_l1 = find((atoms_l1_sym(i1, :) == sa) & (atoms_l1_sym(i2, :) == sa));
        %         else
        %             maxval_l1=-inf;
        %         end
        %         if param.Sfixed
        %             maxval_l1=-inf;
        %         end
        
        
        %%
        if debug
            fprintf('   maxval_l1=%f < %f     maxval_om=%f < %f\n',maxval_l1,param.mu,maxval_om,param.lambda*cf);
            fprintf('   maxval_l1/mu = %f maxval_om/(cf*lambda) = %f\n',maxval_l1/param.mu,maxval_om/(cf*param.lambda));
        end
        
        %% no atoms added, break
        if maxval_l1==-inf && maxval_om==-inf
            fprintf('\n maxval_l1 maxval_om are -inf\n');
            %             keyboard;
            break
        end
        
        %% adding new atom
        %         if (new_row==new_col && maxval_l1/param.mu<=maxval_om/(cf*param.lambda)) || (new_row~=new_col && 2*maxval_l1/param.mu<=maxval_om/(cf*param.lambda))
        %% ading omega atom
        if maxval_om<cf*param.lambda
            fprintf('\n not good atom omega d=%f\n',maxval_om);
            keyboard;
        end
        if debug
            fprintf('\n adding omega atom\n');
            %                 keyboard;
        end
        
        if ActiveSet.atom_count>0
            aom=[ActiveSet.atoms(:,1:ActiveSet.atom_count) anew];
        else
            aom=anew;
        end
        
        %             else
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
                    temp=sanew'*ActiveSet.atoms(:,i);
                    Hall_new(nat+1,i)=temp;
                    Hall_new(i,nat+1)=temp;
                end                
            end
            Hall_new(nat+1,nat+1)=sanew'*anew;
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
        %         else
        %             %% adding l1 atom
        %             if maxval_l1<param.mu
        %                 fprintf('\n not good atom l1 \n');
        %             end
        %             if debug
        %                 fprintf('\n adding l1 atom row=%d col=%d\n', new_row, new_col);
        %                 %                 keyboard;
        %             end
        %             if ~isempty(ActiveSet.I_l1) && sum(ActiveSet.I_l1==idx_l1)
        %                 new_atom_added=false;
        %                 fprintf('\n this l1 atom is already in the collection\n');
        %                 %                 keyboard;
        %                 break;
        %             else
        %                 ActiveSet.I_l1=[ActiveSet.I_l1 idx_l1]; %to avoid adding same atom
        %                 if ActiveSet.atom_count>0
        %                     aom=ActiveSet.atoms(:,1:ActiveSet.atom_count);
        %                 else
        %                     aom=[];
        %                 end
        %
        %                 if debug_update
        %                     if param.f==4
        %                         [Hall_new0,fall_new0] = build_Hessian_l1_sym(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
        %                         [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
        %                     elseif param.f==5
        %                         [Hall_new0,fall_new0] = build_Hessian_l1_SM(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
        %                         [Hall_new,fall_new] = update_Hessian_l1_SM(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
        %                     end
        %
        %                     if norm(Hall_new-Hall_new0,'fro')^2 >1e-10 || norm(fall_new-fall_new0,'fro')^2 >1e-10
        %                         figure(10);clf;
        %                         subplot(1,3,1);
        %                         imagesc(Hall);
        %                         pbaspect([1 1 1]);
        %                         subplot(1,3,2);
        %                         imagesc(Hall_new0);
        %                         pbaspect([1 1 1]);
        %                         subplot(1,3,3);
        %                         imagesc(Hall_new);
        %                         pbaspect([1 1 1]);
        %                         error('the update is not correct\n');
        %                     end
        %                 end
        %                 %                 else
        %                 if param.f==1
        %                     [Hall_new,fall_new] = update_Hessian_prox(inputData.Y,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
        %                 elseif param.f==4
        %                     %                         keyboard;
        %                     [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
        %                 elseif param.f==5
        %                     [Hall_new,fall_new] = update_Hessian_l1_SM(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
        %                 end
        %                 %                 end
        %                 %                 keyboard;
        %
        %                 g=Hall_new*[ActiveSet.beta;0;ActiveSet.alpha]+fall_new;
        %                 idx=length(ActiveSet.beta)+1;
        %                 if g(idx)>0
        %                     fprintf('\n Not a descent direction when adding l1 atom d=%f\n',g(idx));
        %                     ActiveSet.I_l1(end)=[];
        %                     break;
        %                 else
        %                     new_atom_added=true;
        %                     new_atom_l1=true;
        %                     nbetas=nbetas+1;
        %                     Hall=Hall_new;
        %                     fall=fall_new;
        %                 end
        %             end
        %         end
        
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


if 1 %param.debug
    figure(30);clf;
    plot(hist.obj);
    title('objective');
    keyboard;
end

end


