function [ Z,Z1,Z2,Hall,fall, ActiveSet, hist] = solve_ps_l1_omega_asqp( Z,Z1,Z2,ActiveSet,param,inputData,atoms_l1_sym,Hall,fall)
%Using Active Set to solve (PS) problem


fus=false; %fusionning correlated atoms

debug_update=0;
debug=0;
compute_dg=1;
display('in solve_ps change S when changing loss fun\n');
if param.f==4
    S=inputData.X1*inputData.X1;
elseif param.f==5
    S=inputData.X;
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

if param.debug
    alphaSparsity=[];
    diffalphas=[];
end

new_atom_added=false;
idx_added=-1;
first_pass=1;
fusioned_atoms=false;
i=1;
count=1;
cont=true;

H = gradient(Z,inputData,param);
maxIJ=max(abs(H(:)));
maxvar=inf;
res_corr=0; % remove too correlated atoms

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
    
    if ~first_pass %|| ActiveSet.atom_count==0
        
        if new_atom_added
            if new_atom_om
                alpha0=[ActiveSet.beta;ActiveSet.alpha;res_corr];
                idx_added=length(alpha0);
            elseif new_atom_l1
                alpha0=[ActiveSet.beta;0;ActiveSet.alpha];
                idx_added=length(ActiveSet.beta)+1;
            end
        else
            alpha0=[ActiveSet.beta;ActiveSet.alpha];
            idx_added=1;
        end
        
        % active-set
        if length(alpha0)~=size(Hall,1)
            fprintf('sizes mismatch\n');
            keyboard;
        end
        %rcond(Hall+1e-10*eye(length(fall)))
        if debug
            obj0=.5*alpha0'*Hall*alpha0+fall'*alpha0;
        end
        %         fprintf('    nb l1 atoms=%d    nb om atoms=%d\n',length(ActiveSet.I_l1),ActiveSet.atom_count);
        if ActiveSet.atom_count>param.max_nb_atoms
            fprintf('    nb l1 atoms=%d    nb om atoms=%d\n',length(ActiveSet.I_l1),ActiveSet.atom_count);
            fprintf('    max nb om atoms %d reached \n',param.max_nb_atoms);
            cont=false;
            break;
        end
        
        [alph,Jset,npiv]=asqp2(Hall+0*1e-12*eye(length(fall)),-fall,alpha0,param_as,new_atom_added,idx_added);
        %         fprintf('out asqp norm(grad)=%f\n',ng);
        
        
        fusioned_atoms=false;
        
        if debug
            obj1=.5*alph'*Hall*alph+fall'*alph;
            if obj1>obj0
                fprintf('objective increasing in asqp\n');
            end
        end
        
        %% for robustness (too small alpha are deleted)
        eps_alph=1e-8;
        if(min(abs(alph(Jset))))<eps_alph
            fprintf('small alph\n');
            alph(abs(alph)<eps_alph | alph<0)=0;
            Jset(abs(alph)<eps_alph | alph<0)=0;
        end
        
        %%
        if nbetas>0
            Jbeta=Jset(1:nbetas);
            ActiveSet.beta=alph(1:nbetas);
            ActiveSet.beta=ActiveSet.beta(Jbeta);
            ActiveSet.I_l1=ActiveSet.I_l1(Jbeta);
        end
        
        if length(alph)>nbetas
            Jalpha=Jset((nbetas+1):end);
            ActiveSet.alpha=alph((nbetas+1):end);
            %% for robustness, fuisioning too correlated atoms
            if fus
                if ActiveSet.atom_count>1
                    atom=ActiveSet.atoms(:,ActiveSet.atom_count);
                    correl = 1-abs(sum(bsxfun(@times,ActiveSet.atoms(:,1:ActiveSet.atom_count),atom),1));
                    K0= correl'<1e-8;
                    K = Jalpha & K0;
                    if sum(K)>1
                        fprintf('\n too correlated atoms\n');
                        fusioned_atoms=true;
                        %                     keyboard;
                        Jalpha = Jalpha & ~K0;
                        idx=find(K);
                        idx=idx(1);
                        Jalpha(idx)=1;
                        v=sum(full(bsxfun(@times,ActiveSet.atoms(:,K),ActiveSet.alpha(K)')),2);
                        %                     keyboard;
                        ActiveSet.atoms(:,idx)=v/norm(v);
                        ActiveSet.alpha(idx)=norm(v);
                        Jset((nbetas+1):end)=Jalpha;
                    end
                end
            end
            new_atom_count=sum(Jalpha);
            ActiveSet.atom_count=new_atom_count;
            ActiveSet.atoms=ActiveSet.atoms(:,Jalpha);%not necessary (for debbuggging here)
            ActiveSet.alpha=ActiveSet.alpha(Jalpha);
        end
        
        nbetas=length(ActiveSet.beta);
        Hall=Hall(Jset,Jset);
        fall=fall(Jset);
        if debug
            ab=[ActiveSet.beta;ActiveSet.alpha];
            obj2=.5*ab'*Hall*ab+fall'*ab;
        end
        
        %%
        if fusioned_atoms && fus
            %             keyboard;
            alpha0=[ActiveSet.beta;ActiveSet.alpha];
            [alph2,Jset2,npiv]=asqp2(Hall+0*1e-12*eye(length(fall)),-fall,alpha0,param_as,false,idx_added);
            %             keyboard;
            if nbetas>0
                Jbeta=Jset2(1:nbetas);
                ActiveSet.beta=alph2(1:nbetas);
                ActiveSet.beta=ActiveSet.beta(Jbeta);
                ActiveSet.I_l1=ActiveSet.I_l1(Jbeta);
            end
            if length(alph)>nbetas
                ActiveSet.alpha=alph2((nbetas+1):end);
            end
        end
        %%
        
        %% Update ActiveSet and Z
        if 1 %debug
            Z1_old=Z1;
            Z2_old=Z2;
        end
        
        if ~param.Sfixed
            Z1=zeros(p);
            nz=find(ActiveSet.beta>1e-15);
            for j=nz'
                Z1=Z1+ActiveSet.beta(j)*reshape(atoms_l1_sym(:,ActiveSet.I_l1(j)),p,p);
            end
        end
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
                    if obj(i)>obj(i-1)+1e-10
                        fprintf('objective increasing\n');
                        keyboard;
                    end
                end
            end
            i=i+1;
            %% Verify sttopping criterion
            H = gradient(Z,inputData,param);
            [maxIJ,new_row, new_col] = dual_l1_spca(H);
            cf=1;
            if ~isempty(ActiveSet.I)
                [new_i,new_val,maxvar]=get_new_atom_spca(H,ActiveSet,param);
                cf=min(param.cardfun(length(new_i):end));
            else
                maxvar=0;
            end
            %% extra condition
            omega=pen(i-1);
            dotHZ=trace(-H*Z);
            dualomega=max(maxvar/(cf*param.lambda),maxIJ/param.mu);
            cond=omega*dualomega - dotHZ;
            epscond=1e-6;
            %sanity check (output of active set)
            if ~isempty(ActiveSet.I)
                valmax_om=-inf;
                atmax_om=-1;
                for at=1:ActiveSet.atom_count
                    atom=ActiveSet.atoms(:,at);
                    supp=sum(abs(atom)>0);
                    cf=min(param.cardfun(supp:end));
                    val= abs(-atom'*H*atom/cf-param.lambda);
                    if val>valmax_om
                        valmax_om=val;
                        atmax_om=at;
                    end
                end
            else
                valmax_om=0;
            end
            valmax_l1=-inf;
            for at=length(ActiveSet.I_l1)
                atom=atoms_l1_sym(:,ActiveSet.I_l1(at));
                val=-dot(H(:),atom);
                %                 if abs(sum(aom))>1
                %                     val=val/2;
                %                 end
                val=abs(val-param.mu);
                if val>valmax_l1
                    valmax_l1=val;
                end
            end
            if valmax_l1>param.epsStop/2 || valmax_om>param.epsStop/2
                fprintf('Sanity check :  |maxIJ-mu|=%f<%f  |maxvart-lambda|=%f<%f \n',valmax_l1,param.epsStop,valmax_om,param.epsStop);
                %                 keyboard;
            end
            % end of sanity check (output of active set
            
            % Stopping criterion
%             if maxvar < param.lambda*(1+param.epsStop) && maxIJ < param.mu*(1+param.epsStop) && cond<epscond %&& dualgap/param.PSdualityEpsilon<1
%                 cont=false;
%             end
            if cond<epscond %&& dualgap/param.PSdualityEpsilon<1
                cont=false;
            end
            fprintf('maxIJ/mu=%4.2f<1     varmax/cf*lambda=%4.2f<1   dg/eps=%4.2f<1  cond=%4.2f<%4.2f\n',maxIJ/param.mu,maxvar/(cf*param.lambda),dualgap/param.PSdualityEpsilon,cond,epscond);
            if debug
                fprintf('  maxIJ/mu=%4.2f < 1     varmax/lambda=%4.2f < 1 var-varold=%4.2f continue=%d  count=%d\n',maxIJ/param.mu,maxvar/param.lambda,norm(maxvar-maxvarold),cont && count< param.niterPS,count);
                %                     keyboard;
            end
            
            cont=cont && count< param.niterPS;
        end
    end
    
    %% get new atom
    if debug
        fprintf('\n--------------------------------------------\n');
    end
    
    if cont
        
        new_atom_l1=false;
        new_atom_om=false;
        H = gradient(Z,inputData,param);
        
        %% omega atom
        
        cf=1;
        if ~isempty(ActiveSet.I)
            [new_i, new_val, maxval_om0]=get_new_atom_spca(H,ActiveSet,param);
            anew=sparse(new_i,ones(length(new_i),1),new_val,p,1);
            cf=min(param.cardfun(length(new_i):end));
%             if maxval_om0>param.lambda*(1+param.epsStop)
            if maxval_om0-cf*param.lambda>1e-16
                maxval_om=maxval_om0;
            else
                maxval_om=-inf; 
            end
        else
            maxval_om=-inf;
        end
        
        %% l1 sym atom
        
        [maxval_l1,new_row, new_col] = dual_l1_spca(H);
        
%         if maxval_l1>param.mu*(1+param.epsStop)
        if maxval_l1-param.mu>1e-16
            sa=-sign(H(new_row,new_col));
            i1=(new_row-1)*p+new_col;
            i2=(new_col-1)*p+new_row;
            idx_l1 = find((atoms_l1_sym(i1, :) == sa) & (atoms_l1_sym(i2, :) == sa));
        else
            maxval_l1=-inf;
        end
        if param.Sfixed
            maxval_l1=-inf;
        end
        
        
        %%
        if debug
            fprintf('  maxval_l1=%f < %f     maxval_om=%f < %f\n',maxval_l1,param.mu,maxval_om,param.lambda*cf);
            fprintf('maxval_l1/mu = %f maxval_om/lambda = %f\n',maxval_l1/param.mu,maxval_om/param.lambda);
        end
        
        %% no atoms added, break
        if maxval_l1==-inf && maxval_om==-inf
            fprintf('\n maxval_l1 maxval_om are -inf\n');
            %             keyboard;
            break
        end
        
        %% adding new atom
        if (new_row==new_col && maxval_l1/param.mu<=maxval_om/(cf*param.lambda)) || (new_row~=new_col && 2*maxval_l1/param.mu<=maxval_om/(cf*param.lambda))
            %% ading omega atom
            if maxval_om<cf*param.lambda
                fprintf('\n not good atom omega d=%f\n',maxval_om);
                keyboard;
            end
            if debug
                fprintf('\n adding omega atom\n');
            end
            
            if ActiveSet.atom_count>0
                aom=[ActiveSet.atoms(:,1:ActiveSet.atom_count) anew];
            else
                aom=anew;
            end
            if debug_update
                if param.f==4
                    [Hall_new0,fall_new0] = build_Hessian_l1_sym(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
                    [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
                elseif param.f==5
                    [Hall_new0,fall_new0] = build_Hessian_l1_SM(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
                    [Hall_new,fall_new] = update_Hessian_l1_SM(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
                end
                if norm(Hall_new-Hall_new0,'fro')^2 >1e-10 || norm(fall_new-fall_new0,'fro')^2 >1e-10
                    figure(10);clf;
                    subplot(1,3,1);
                    imagesc(Hall);
                    pbaspect([1 1 1]);
                    subplot(1,3,2);
                    imagesc(Hall_new0);
                    pbaspect([1 1 1]);
                    subplot(1,3,3);
                    imagesc(Hall_new);
                    pbaspect([1 1 1]);
                    error('the update is not correct\n');
                end
            else
                if param.f==4
                    [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
                elseif param.f==5
                    [Hall_new,fall_new] = update_Hessian_l1_SM(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
                end
            end
            
            g=Hall_new*[ActiveSet.beta;ActiveSet.alpha;0]+fall_new;
            if g(end)>0
                fprintf('\n Not a descent direction when adding om atom d=%f\n',g(end));
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
        else
            %% adding l1 atom
            if maxval_l1<param.mu
                fprintf('\n not good atom l1 \n');
            end
            if debug
                fprintf('\n adding l1 atom row=%d col=%d\n', new_row, new_col);
            end
            if ~isempty(ActiveSet.I_l1) && sum(ActiveSet.I_l1==idx_l1)
                new_atom_added=false;
                fprintf('\n this l1 atom is already in the collection\n');
                %                 keyboard;
                break;
            else
                ActiveSet.I_l1=[ActiveSet.I_l1 idx_l1]; %to avoid adding same atom
                if ActiveSet.atom_count>0
                    aom=ActiveSet.atoms(:,1:ActiveSet.atom_count);
                else
                    aom=[];
                end
                
                if debug_update
                    if param.f==4
                        [Hall_new0,fall_new0] = build_Hessian_l1_sym(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
                        [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
                    elseif param.f==5
                        [Hall_new0,fall_new0] = build_Hessian_l1_SM(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
                        [Hall_new,fall_new] = update_Hessian_l1_SM(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
                    end
                    
                    if norm(Hall_new-Hall_new0,'fro')^2 >1e-10 || norm(fall_new-fall_new0,'fro')^2 >1e-10
                        figure(10);clf;
                        subplot(1,3,1);
                        imagesc(Hall);
                        pbaspect([1 1 1]);
                        subplot(1,3,2);
                        imagesc(Hall_new0);
                        pbaspect([1 1 1]);
                        subplot(1,3,3);
                        imagesc(Hall_new);
                        pbaspect([1 1 1]);
                        error('the update is not correct\n');
                    end
                else
                    if param.f==4
                        [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
                    elseif param.f==5
                        [Hall_new,fall_new] = update_Hessian_l1_SM(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
                    end
                end
                %                 keyboard;
                
                g=Hall_new*[ActiveSet.beta;0;ActiveSet.alpha]+fall_new;
                idx=length(ActiveSet.beta)+1;
                if g(idx)>0
                    fprintf('\n Not a descent direction when adding l1 atom d=%f\n',g(idx));
                    ActiveSet.I_l1(end)=[];
                    break;
                else
                    new_atom_added=true;
                    new_atom_l1=true;
                    nbetas=nbetas+1;
                    Hall=Hall_new;
                    fall=fall_new;
                end
            end
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
    figure(30);clf;
    subplot(1,2,1);
    plot(hist.obj);
    title('objective');
    subplot(1,2,2);
    semilogy(hist.dg);
    title('duality gap');
    figure(31);clf
    subplot(1,3,1);
    imagesc(abs(Z1))
    pbaspect([1 1 1]);
    subplot(1,3,2);
    imagesc(abs(Z2))
    pbaspect([1 1 1]);
    subplot(1,3,3);
    imagesc(abs(Z))
    pbaspect([1 1 1]);
    keyboard;
end

end


