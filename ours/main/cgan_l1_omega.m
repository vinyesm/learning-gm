function [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param,startingZ,ActiveSet)

MAX_NB_ATOMS=50;
param.max_nb_atoms=MAX_NB_ATOMS;

%% init
if nargin < 3
    startingZ = set_default_Z(inputData,param);
    Z = startingZ;
    Z1 =Z;
    Z2=Z;
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
else
    Z = startingZ;
end

param = set_default_param(param);

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

fprintf('Warning : change build_atoms_hessian_l1_sym when loss is not .5*|S^.5*X*S.^5-I|\n');
[ Q,q,atoms_l1_sym ] = build_atoms_hessian_l1_sym(inputData.X1*inputData.X1,param.mu);
% ActiveSet.beta=zeros(size(Q,1),1);
% ActiveSet.I_l1=1:size(Q,1);
% Hall=Q;
% fall=q+param.mu*ones(size(Q,1),1);
ActiveSet.beta=[];
ActiveSet.I_l1=[];
Hall=[];
fall=[];
cardVal=[];
U=[];

tic

max_nb_atoms = param.max_nb_atoms;
max_nb_main_loop = param.max_nb_main_loop;
hist.time=toc;

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
        
        [Z, Z1, Z2,U,Hall,fall,cardVal, ActiveSet, hist_ps] = solve_ps_l1_omega_asqp(Z,Z1,Z2, ActiveSet,param,inputData,atoms_l1_sym,U,Hall,fall,cardVal);
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
    
    %% get a new descent direction using truncated power iteration
    
    H = gradient(Z,inputData,param);
%     keyboard;
    
    if param.verbose==1
        fprintf('%d/%d   \n',i,max_nb_main_loop);
    end
    
    [u, kBest,val] = lmo_spsd_TPower(-H,param);
    if val<0
        currI=[];
        fprintf('   all eigs are negative\n')
%         keyboard;
    else
        param.k=kBest;
        currI = find(u);
    end
%     keyboard;
    
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
        fprintf('   varIJ = %2.4e, thresh = %2.4e, length(currI)=%d\n',varIJ, param.lambda*(1+param.epsStop / kBest)* param.cardfun(kBest), length(currI))
    end
    
    
    if varIJ < param.lambda*(1+param.epsStop) && maxIJ < param.mu*(1+param.epsStop)
        c=0;
    elseif ActiveSet.atom_count>=param.max_nb_atoms
%         keyboard;
        c=0;
    elseif takenI
        fprintf('This support has already been added. Stopping\n');
        %c=0;
    elseif varIJ > param.lambda*(1+param.epsStop)
        ActiveSet.I = [ActiveSet.I, currI];
        %ActiveSet.U = [ActiveSet.U, u(currI)];
        ActiveSet.Sigma = [ActiveSet.Sigma, varIJ];
        ActiveSet.Z = [ActiveSet.Z, zeros(param.k,param.k)];
        ActiveSet.tracenorm = [ ActiveSet.tracenorm , 0];
        ActiveSet.k = [ActiveSet.k , kBest];
        ActiveSet.fronorm = [ ActiveSet.fronorm , 0];
    else
        %c = 0;
    end
    c = i<max_nb_main_loop & c;
end

if param.debug==1
    if i>=max_nb_main_loop
        fprintf('\n max number of main loop iterations reached\n');
    end
end

% % ActiveSet = postProcessFactors(ActiveSet,Z);
% 
% al1=atoms_l1_sym(:,ActiveSet.I_l1);
% aom=ActiveSet.atoms;
% [Hp,fp] = build_Hessian_postprocess(inputData,param,al1,aom);
% lb=zeros(1,length(fp));
% ub=ones(1,length(fp));
% param_as.max_iter=1e3;
% param_as.epsilon=1e-14;
% param_as.debug_mode=false;
% param_as.ws=true;
% 
% %options.OptimalityTolerance=1e-10;
% [alph]=quadprog(Hp,fp,[],[],[],[],lb,ub);
% keyboard;
% Jset= alph>1e-2;
% 
% nbetas=length(ActiveSet.beta);
% if length(alph)>nbetas
%     Jalpha=Jset((nbetas+1):end);
%     new_atom_count=sum(Jalpha);
%     ActiveSet.atom_count=new_atom_count;
%     ActiveSet.atoms=ActiveSet.atoms(:,Jalpha);%not necessary (for debbuggging here)
%     ActiveSet.alpha=alph((nbetas+1):end);
%     ActiveSet.alpha=ActiveSet.alpha(Jalpha);
%     %             U=U(:,Jalpha);
%     cardVal=cardVal(Jalpha);
% end
% %ActiveSet.alpha
        


%%
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

end






