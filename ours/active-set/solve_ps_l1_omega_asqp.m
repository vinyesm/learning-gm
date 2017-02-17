function [ Z,Z1,Z2,U,Hall,fall,cardVal, ActiveSet, hist] = solve_ps_l1_omega_asqp( Z,Z1,Z2,ActiveSet,param,inputData,atoms_l1_sym,U,Hall,fall,cardVal)
%Using Active Set to solve (PS) problem

debug_update=0;
debug=0;
compute_dg=0;
display('in solve_ps change S when changing loss fun\n');
S=inputData.X1*inputData.X1;
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
i=1;
count=1;
cont=true;

H = gradient(Z,inputData,param);
maxII=max(abs(H(:)));
maxvar=inf;

while cont
    %% get a new atom ui  [si=ui'*ui] in ActiveSet.atoms  from some C_I I
    % in ActiveSet.I to add to the collection
    % current point is Z
    
    if ActiveSet.atom_count>param.max_nb_atoms
        display('maximum number of atoms added');
        break;
    end
    
    if ~first_pass %|| ActiveSet.atom_count==0
        
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
        
        % active-set
        if length(alpha0)~=size(Hall,1)
            fprintf('sizes mismatch\n');
            keyboard;
        end
        %rcond(Hall+1e-10*eye(length(fall)))
        if debug
            obj0=.5*alpha0'*Hall*alpha0+fall'*alpha0;
        end
        [alph,Jset,npiv]=asqp2(Hall+0*1e-12*eye(length(fall)),-fall,alpha0,param_as,new_atom_added,idx_added);
        if debug
            obj1=.5*alph'*Hall*alph+fall'*alph;
            if obj1>obj0
                fprintf('objective increasing in asqp\n');
                %                 keyboard;
            end
        end
        
        eps_alph=1e-10;
        if(min(abs(alph(Jset))))<eps_alph
            fprintf('small alph\n');
            alph(abs(alph)<eps_alph || alph<0)=0;
            Jset(abs(alph)<eps_alph || alph<0)=0;
            %             keyboard;
        end
        
        if nbetas>0
            Jbeta=Jset(1:nbetas);
            ActiveSet.beta=alph(1:nbetas);
            ActiveSet.beta=ActiveSet.beta(Jbeta);
            ActiveSet.I_l1=ActiveSet.I_l1(Jbeta);
            
        end
        
        if length(alph)>nbetas
            Jalpha=Jset((nbetas+1):end);
            new_atom_count=sum(Jalpha);
            ActiveSet.atom_count=new_atom_count;
            ActiveSet.atoms=ActiveSet.atoms(:,Jalpha);%not necessary (for debbuggging here)
            ActiveSet.alpha=alph((nbetas+1):end);
            ActiveSet.alpha=ActiveSet.alpha(Jalpha);
            %             U=U(:,Jalpha);
            cardVal=cardVal(Jalpha);
        end
        %ActiveSet.alpha
        nbetas=length(ActiveSet.beta);
        Hall=Hall(Jset,Jset);
        fall=fall(Jset);
        if debug
            ab=[ActiveSet.beta;ActiveSet.alpha];
            obj2=.5*ab'*Hall*ab+fall'*ab;
        end
        
        
        %% Update ActiveSet and Z
        Z1=zeros(p);
        nz=find(ActiveSet.beta>1e-15);
        for j=nz'
            Z1=Z1+ActiveSet.beta(j)*reshape(atoms_l1_sym(:,ActiveSet.I_l1(j)),p,p);
        end
        Z2=zeros(p);
        nz=find(ActiveSet.alpha>1e-15);
        for j=nz'
            u=ActiveSet.atoms(:,j);
            Z2=Z2+ActiveSet.alpha(j)*(u*u');
        end
        Z=Z1+Z2;
        
        if debug
            obj3=.5*norm(inputData.Y - inputData.X1*Z*inputData.X2, 'fro')^2+param.lambda*sum(ActiveSet.alpha)+param.mu*sum(ActiveSet.beta)-.5*p;
            if obj3-obj2>1e-10
                fprintf('objective increasing after constructing Z \n');
                keyboard;
            end
        end
        
        
        %% Compute objective, loss, penalty and duality gap
        if (param.sloppy==0 || (param.sloppy~=0 && mod(count,100)==1)) %&& ~isempty(ActiveSet.alpha)
            if compute_dg
            [loss(i),pen(i),obj(i),dg(i),time(i)]=get_val_l1_omega_asqp(Z,ActiveSet,inputData,param,cardVal);
            nb_pivot(i)=npiv;
            active_var(i)= sum(ActiveSet.alpha>0);
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
            maxII=max(abs(H(:)));
            if ~isempty(ActiveSet.I)
                maxvarold=maxvar;
                [maxvar, kmaxvar]=max_var(Z,ActiveSet,param,inputData );
                if debug && norm(maxvar-maxvarold)<1e-10
                    fprintf('maxvar not changing\n');
                    keyboard;
                end
                size_supp=length(ActiveSet.I{kmaxvar});
                if maxvar < param.lambda*(1+param.epsStop / size_supp)* param.cardfun(size_supp) && maxII<param.mu*(1+param.epsStop)
                    cont=false;
                end
                if debug 
                    fprintf('  maxII=%f < %f     varmax=%f < %f\n',maxII,param.mu,maxvar,param.lambda);
                end
            else
                if  maxII<param.mu*(1+param.epsStop)
                    cont=false;
                end
                if debug 
                    fprintf('maxII=%f < %f \n',maxII,param.mu);
                end
            end
            %             if count> param.niterPS;
            %                 keyboard;
            %             end
            cont=cont && count< param.niterPS;
        end
    end
    
    %     if isempty(ActiveSet.I)
    %         cont=false;
    %     end
    
    %% get new atom
    if debug
        fprintf('\n--------------------------------------------\n');
    end
    if cont
        
        new_atom_l1=false;
        new_atom_om=false;
        
        %% omega atom
        
        if ~isempty(ActiveSet.I)
            StartY=inputData.Y;
            inputData.Y=StartY-inputData.X1*Z1*inputData.X2;
            [new_i, new_val, maxval_om]=get_new_atom_spca(Z2,ActiveSet,param,inputData);
            inputData.Y=StartY;
            if maxval_om<param.lambda*(1+param.epsStop)
                maxval_om=-inf;
            end
        else
            maxval_om=-inf;
        end
        
        %% l1 sym atom
        
        if maxII>param.mu*(1+param.epsStop)
            StartY=inputData.Y;
            inputData.Y=StartY-inputData.X1*Z2*inputData.X2;
            H = gradient(Z1,inputData,param);
            [maxval_l1]=max(abs(H(:)));
            %         [maxval_l1]=max((-H(:)));
            inputData.Y=StartY;
            
            [new_row, new_col] = find(abs(H) == maxval_l1);
            new_row=new_row(1);
            new_col=new_col(1);
            sa=-sign(maxval_l1*H(new_row,new_col));
            i1=(new_row-1)*p+new_col;
            i2=(new_col-1)*p+new_row;
            idx_l1 = find((atoms_l1_sym(i1, :) == sa) & (atoms_l1_sym(i2, :) == sa));
        else
            maxval_l1=-inf;
        end
        
        %%
        if debug
            fprintf('  maxII=%f < %f     maxval_om=%f < %f\n',maxII,param.mu,maxval_om,param.lambda);
            fprintf('maxval_l1/mu = %f maxval_om/lambda = %f\n',maxval_l1/param.mu,maxval_om/param.lambda);
%             if maxval_om/param.lambda~=-inf && maxval_om/param.lambda<0
%                 keyboard
%             end
        end
        if maxval_l1==-inf && maxval_om==-inf
%             keyboard;
            break
        end
        
        if maxval_l1/param.mu<=maxval_om/param.lambda
            %ading omega atom
            if maxval_om<param.lambda
                fprintf('\n not good atom omega d=%f\n',maxval_om);
                keyboard;
            end
            if debug
                fprintf('\n adding omega atom\n');
            end
            anew=sparse(new_i,ones(length(new_i),1),new_val,p,1);
            %[Hall_new,fall_new,U_new] = add_omega_atoms_hessian_l1_sym(inputData,param,Hall,fall,atoms_l1_sym,U,anew);
            if ActiveSet.atom_count>0
                aom=[ActiveSet.atoms(:,1:ActiveSet.atom_count) anew];
            else
                aom=anew;
            end
            if debug_update
                [Hall_new0,fall_new0] = build_Hessian_l1_sym(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
                [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
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
                [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,2);
            end
            
            g=Hall_new*[ActiveSet.beta;ActiveSet.alpha;0]+fall_new;
            if g(end)>0
                fprintf('\n Not a descent direction d=%f\n',g(end));
                %                 keyboard;
                break;
            else
                ActiveSet.atom_count = ActiveSet.atom_count +1;
                ActiveSet.max_atom_count_reached=max(ActiveSet.max_atom_count_reached,ActiveSet.atom_count);
                ActiveSet.atomsSupport=[ActiveSet.atomsSupport new_i];
                ActiveSet.atoms(:,ActiveSet.atom_count)=anew;
                Hall=Hall_new;
                fall=fall_new;
                %U=U_new;
                weight=1;
                cardVal=[cardVal;weight];
                new_atom_added=true;
                new_atom_om=true;
            end
        else
            %adding l1 atom
            if maxval_l1<param.mu
                fprintf('\n not good atom l1 \n',maxval_l1);
                %                 keyboard;
            end
            if debug
                fprintf('\n adding l1 atom row=%d col=%d\n', new_row, new_col);
            end
            if sum(ActiveSet.I_l1==idx_l1)
                new_atom_added=false;
                fprintf('\n this atom is already in the collection\n');
                %                 keyboard;
            else
                ActiveSet.I_l1=[ActiveSet.I_l1 idx_l1]; %to avoid adding same atom
                if ActiveSet.atom_count>0
                    aom=ActiveSet.atoms(:,1:ActiveSet.atom_count);
                else
                    aom=[];
                end
                
                if debug_update
                    [Hall_new0,fall_new0] = build_Hessian_l1_sym(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),aom);
                    [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
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
                    [Hall_new,fall_new] = update_Hessian_l1_sym(S,param,Hall, fall,atoms_l1_sym(:,ActiveSet.I_l1),aom,1);
                end
                
                g=Hall_new*[ActiveSet.beta;0;ActiveSet.alpha]+fall_new;
                idx=length(ActiveSet.beta)+1;
                if g(idx)>0
                    fprintf('\n Not a descent direction d=%f\n',g(idx));
                    ActiveSet.I_l1(end)=[];
                    %                     keyboard;
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


