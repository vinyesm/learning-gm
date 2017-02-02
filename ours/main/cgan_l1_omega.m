function [Z ActiveSet hist param flag output] = cgan_l1_omega(inputData,param,startingZ,ActiveSet)

%% init
if nargin < 3
    startingZ = set_default_Z(inputData,param);
    Z = startingZ;
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
    
    if ~isempty(ActiveSet.I)
        if param.verbose==1
            fprintf('    solving PS..\n ')
        end
        if strcmp('asqp',param.opt)
            [Z, ActiveSet, hist_ps] = solve_ps_l1_omega_asqp(Z, ActiveSet,param,inputData);
            if ~isempty(ActiveSet.alpha) && param.debug==1
                hist.nzalphas=[hist.nzalphas full(sum(ActiveSet.alpha>0))];
                hist.normalpha=[hist.normalpha full(sum(ActiveSet.alpha))];
                fprintf('   nz alphas=%d  |alpha|=%f  dual gap =%f\n',full(sum(ActiveSet.alpha>0)),full(sum(ActiveSet.alpha)), dg(end));
            end
            nb_pivot=[nb_pivot hist_ps.nb_pivot];
            active_var=[active_var hist_ps.active_var];
        end
        if strcmp('proxbcd',param.opt)
            [Z, D, ActiveSet, hist_ps] = solve_ps_spca_proxbcd(Z,D,ActiveSet,param,inputData);
        end
        
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

    if param.verbose==1
        fprintf('%d/%d   \n',i,max_nb_main_loop);
    end
    
    [u, kBest] = lmo_spsd_TPower(-H,param);
    param.k=kBest;
    currI = find(u);

    %% verbose
    if param.verbose==1
        fprintf('   currI = ')
        for j=1:length(currI)
            fprintf('%d ',currI(j));
        end
        fprintf('\n');
        
        if(isempty(currI))
            fprintf('currI is empty\n');
        end
    end
    %%
    
    varIJ = norm(H(currI,currI));
    takenI= isInCell(currI,ActiveSet.I,cell2mat(ActiveSet.k)) ;
    hist.varIJ=[hist.varIJ varIJ];
    
    flag.var=varIJ;
    
    if param.verbose==1
        fprintf('   variance = %2.4e, thresh = %2.4e, length(currI)=%d\n',varIJ, param.lambda*(1+param.epsStop / kBest)* param.cardfun(kBest), length(currI))
    end
    
    
    if takenI
        fprintf('This support has already been added. Stopping\n');
        c=0;
    elseif varIJ > param.lambda*(1+param.epsStop / kBest)* param.cardfun(kBest)
        ActiveSet.I = [ActiveSet.I, currI];
        %ActiveSet.U = [ActiveSet.U, u(currI)];
        ActiveSet.Sigma = [ActiveSet.Sigma, varIJ];
        ActiveSet.Z = [ActiveSet.Z, zeros(param.k,param.k)];
        ActiveSet.tracenorm = [ ActiveSet.tracenorm , 0];
        ActiveSet.k = [ActiveSet.k , kBest];
        ActiveSet.fronorm = [ ActiveSet.fronorm , 0];
    else
        c = 0;
    end
    c = i<max_nb_main_loop & c;
end

if param.debug==1
    if i>=max_nb_main_loop
        fprintf('\n max number of main loop iterations reached\n');
    end
end

% ActiveSet = postProcessFactors(ActiveSet,Z);
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





