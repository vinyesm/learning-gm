function [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param,startingZ,ActiveSet)

pm=1; %postprocessing after PS if no new atom found
pp=1; %postprocessing
pt=0; %postprocessing thresh
MAX_NB_ATOMS=500;
param.max_nb_atoms=MAX_NB_ATOMS;


param = set_default_param(param);
%% init
if nargin < 3
    startingZ = set_default_Z(inputData,param);
    Z = startingZ;
    Z1 = Z;
    Z2 = Z;
    ActiveSet = {};
    ActiveSet.I = {};
    ActiveSet.U = {};
    ActiveSet.Sigma = {};
    ActiveSet.Z = {};
    ActiveSet.tracenorm = {};
    ActiveSet.fronorm = {};
    ActiveSet.k = {};
    ActiveSet.atomsSupport = {};
    ActiveSet.alpha= [];
    ActiveSet.atom_count = 0;
    ActiveSet.atoms=[];
    ActiveSet.max_atom_count_reached=0;
    ActiveSet.beta=[];
    ActiveSet.I_l1=[];
    Hall=[];
    fall=[];
    cardVal=[];
    U=[];
    qs=3:-1:0;
else
    %     Z = startingZ;
    Z1 = startingZ.Z1;
    Z2 = startingZ.Z2;
    Z = Z1+Z2;
    qs=0;
end

if param.Sfixed
    Z1 = param.Sstar;
    Z = Z1;
end


obj  = [];
loss = [];
pen  = [];
dg   = [];
time =[];
dg_sup = [];
time_sup =[];
obj_sup = [];
nb_pivot=[];
active_var=[];
hist.varIJ=[];
hist.nzalphas=[];
hist.normalpha=[];
p=size(Z,1);
D=zeros(p,1);
output.time=0;

fprintf('\n \n CGAN L1 OMEGA\n');

if param.f==1
    [ Q,q,atoms_l1_sym ] = build_atoms_hessian_prox(inputData.Y,param.mu);
elseif param.f==4
%     fprintf('Warning : change build_atoms_hessian_l1_sym when loss is not .5*|S^.5*X*S.^5-I|\n');
    [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_sym(inputData.X1*inputData.X1,param.mu);
elseif param.f==5
    [ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_SM(inputData.X,param.mu); %score matching
end

if nargin > 2
    if param.f==4
        [Hall,fall] = build_Hessian_l1_sym(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),ActiveSet.atoms);
    elseif param.f==5
        [Hall,fall] = build_Hessian_l1_SM(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),ActiveSet.atoms);
    end
end

% keyboard;

tic

max_nb_atoms = param.max_nb_atoms;
max_nb_main_loop = param.max_nb_main_loop;
hist.time=toc;
epsStop=param.epsStop;



% for q=5:-1:0
for q=qs
    param.epsStop=2^q*epsStop;
    c = 1;
    i = 0;
    while c
        i = i+1;
        
        if ActiveSet.atom_count > max_nb_atoms
            break;
        end
        
        %% solve problem P_S
        
        if 1 %~isempty(ActiveSet.I) % now we already have all the Eij+Eji atoms
            if param.verbose==1
                fprintf('   solving PS..\n ')
            end
            
            param.epsStop=2^(q-1)*epsStop;
            %             keyboard;
            [Z, Z1, Z2,Hall,fall, ActiveSet, hist_ps] = solve_ps_l1_omega_asqp(Z,Z1,Z2, ActiveSet,param,inputData,atoms_l1_sym,Hall,fall);
            %             keyboard;
            param.epsStop=2^q*epsStop;
            
            if ~isempty(ActiveSet.alpha) && param.debug==1
                hist.nzalphas=[hist.nzalphas full(sum(ActiveSet.alpha>0))];
                hist.normalpha=[hist.normalpha full(sum(ActiveSet.alpha))];
                fprintf('   nz alphas=%d  |alpha|=%f  dual gap =%f\n',full(sum(ActiveSet.alpha>0)),full(sum(ActiveSet.alpha)), hist_ps.dg(end));
            end
            nb_pivot=[nb_pivot hist_ps.nb_pivot];
            active_var=[active_var hist_ps.active_var];
            
            
            obj = [obj hist_ps.obj];
            loss = [loss hist_ps.loss];
            pen = [pen hist_ps.pen];
            dg = [dg hist_ps.dg];
            time = [time hist_ps.time];
            dg_sup = [dg_sup hist_ps.dg_sup];
            time_sup = [time_sup hist_ps.time_sup];
            obj_sup = [obj_sup hist_ps.obj_sup];
        end
        
        %% Cleaning, proximal steps
        
        %% get a new descent direction using truncated power iteration
        
        H = gradient(Z,inputData,param);
        
        if param.verbose==1
            fprintf('%d/%d   \n',i,max_nb_main_loop);
        end
        
        [u, kBest,val] = lmo_spsd_TPower(-H,param);
        cf=min(param.cardfun(kBest:end));
        %     keyboard;
        
        if val<param.lambda*(1+param.epsStop)
            %%      few proximal steps for postprcessing
            %             keyboard;
            if pm && ActiveSet.atom_count>0 && param.f~=1
                fprintf('No new atom found, prox steps for cleaning after PS.. \n');
                if param.f==4
                    S=inputData.X1*inputData.X1;
                elseif param.f==5
                    S=inputData.X;
                end
                [Z2,ActiveSet]=prox_cleaning(Z1,Z2,S,ActiveSet,param,20,0);
                Z=Z1+Z2;
                if param.f==4
                    [Hall,fall] = build_Hessian_l1_sym(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),ActiveSet.atoms);
                elseif param.f==5
                    [Hall,fall] = build_Hessian_l1_SM(inputData,param,atoms_l1_sym(:,ActiveSet.I_l1),ActiveSet.atoms);
                end
                H = gradient(Z,inputData,param);
                val_old=val;
                [u, kBest,val] = lmo_spsd_TPower(-H,param);
                cf=min(param.cardfun(kBest:end));
                %             fprintf('old val=%f new val=%f < %f.. \n',val_old,val, param.lambda*cf);
                %                             keyboard;
            end
        end
        
        if val<0
            currI=[];
            fprintf('   all eigs are negative\n')
            %         keyboard;
        else
            param.k=kBest;
            currI = find(u);
        end
        %         keyboard;
        
        %% verbose
        if param.verbose==1
            fprintf('   lambda=%f  mu=%f\n', param.lambda, param.mu)
            fprintf('   currI = ')
            for j=1:length(currI)
                fprintf('%d ',currI(j));
            end
            fprintf('\n');
            
            if(isempty(currI))
                fprintf('   currI is empty\n');
            end
        end
        %%
        %maxIJ=max(abs(H(:)));
        maxIJ = dual_l1_spca(H);
        if(isempty(currI))
            varIJ=-1;
            takenI=true;
        else
            %         varIJ = norm(H(currI,currI));
            %         varIJ = abs(eigs(H(currI,currI),1,'lm'));
            %fprintf ('TO CHECK: changing stopping criterion to operator norm on currI instead of Frobenius\n')
            varIJ=val;
            takenI= isInCell(currI,ActiveSet.I,cell2mat(ActiveSet.k)) ;
        end
        
        hist.varIJ=[hist.varIJ varIJ];
        flag.var=varIJ;
        
        if param.verbose==1
            fprintf('   maxIJ = %2.4e, thresh = %2.4e\n',maxIJ, param.mu*(1+param.epsStop));
            fprintf('   varIJ = %2.4e, thresh = %2.4e, length(currI)=%d\n',varIJ, param.lambda*(1+param.epsStop / kBest)* param.cardfun(kBest), length(currI));
            rho=p/2;
            if ~isempty(pen)
                dg1=abs(trace(H*Z)+pen(end));
            else
                dg1=abs(trace(H*Z));
            end
            dg2=max(maxIJ-param.mu, varIJ-param.lambda)*rho;
            fprintf('   dg1 = %2.4e dg2 = %2.4e  dg =  %2.4e\n',dg1,dg2,dg1+dg2);
        end
        
        
        if varIJ < param.lambda*cf*(1+param.epsStop) && maxIJ < param.mu*(1+param.epsStop)
            c=0;
        elseif ActiveSet.atom_count>=param.max_nb_atoms
            %         keyboard;
            c=0;
        elseif takenI
            fprintf('This support has already been added. Stopping\n');
            %c=0;
        elseif varIJ > param.lambda*cf*(1+param.epsStop)
            ActiveSet.I = [ActiveSet.I, currI];
            %ActiveSet.U = [ActiveSet.U, u(currI)];
            %ActiveSet.Sigma = [ActiveSet.Sigma, varIJ];
            %ActiveSet.Z = [ActiveSet.Z, zeros(param.k,param.k)];
            %ActiveSet.tracenorm = [ ActiveSet.tracenorm , 0];
            ActiveSet.k = [ActiveSet.k , kBest];
            %ActiveSet.fronorm = [ ActiveSet.fronorm , 0];
        else
            %c = 0;
        end
        c = i<max_nb_main_loop & c;
        %         keyboard
    end
end

if param.debug==1
    if i>=max_nb_main_loop
        fprintf('\n max number of main loop iterations reached\n');
    end
end

%%
if param.Sfixed
    Z1 = param.Sstar;
    Z=Z1+Z2;
end
hist.obj = obj;
hist.loss = loss;
hist.pen = pen;
hist.dg = dg;
hist.time = time;
hist.dg_sup = dg_sup;
hist.time_sup= time_sup;
hist.obj_sup = obj_sup;
hist.nb_pivot= nb_pivot;
hist.active_var= active_var;
if ActiveSet.atom_count>0
    ActiveSet.atoms=ActiveSet.atoms(:,1:ActiveSet.atom_count);
end

%% from atoms to matrices
ActiveSet.matrix_atoms={};
for i=1:ActiveSet.atom_count
    ActiveSet.matrix_atoms{i}=ActiveSet.atoms(:,i)*ActiveSet.atoms(:,i)';
end

%% postprocessing to blocks

if param.f~=1
if pp==1
    if ~isempty(ActiveSet.atoms)
        fprintf('Postprocessing.. \n');
        if param.f==4
            S=inputData.X1*inputData.X1;
        elseif param.f==5
            S=inputData.X;
        end
        [Z2,ActiveSet]=prox_cleaning(Z1,Z2,S,ActiveSet,param,100,1);
        Z=Z1+Z2;
    end
end
end

if pt==1
    if ~isempty(ActiveSet.atoms)
        fprintf('Postprocessing.. \n');
        thresh=1e-6;
        [Z2,ActiveSet]=postprocessing(ActiveSet, thresh);
        Z=Z1+Z2;
    end
end

end






