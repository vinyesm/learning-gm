function [M,S,L,U,hist] = logdetPPA_l1_omega(Sigma,param,ActiveSet)

param.epsStop=1e-8;
param.debug=0;


obj_sup = [];
hist.varIJ =[];
fprintf('\n \n logdetPPA L1 OMEGA\n');

%% set up SDP data in SDPT3 format
%%
p=size(Sigma,1);
k=param.k;


p2 = p*(p+1)/2;
b = zeros(p2,1);
C{1} = Sigma;
blk{1,1} = 's'; blk{1,2} = p;
[Iall,Jall] = find(triu(ones(p)));
tmp = [Iall,Jall];
m2 = size(tmp,1);
Icomp = tmp(:,1); Jcomp = tmp(:,2);
Itmp  = Icomp + Jcomp.*(Jcomp-1)/2;
Atmp  = spconvert([Itmp,[1:m2]',ones(m2,1); p2,m2,0]);
At{1} = Atmp;

blk{2,1} = 'l'; blk{2,2} = 2*p2;
Identity = speye(p2);
At{2,1} = [-Identity,Identity]';
idx = find(Icomp == Jcomp);
ee  = sqrt(2)*ones(m2,1);
if ~isempty(idx); ee(idx) = ones(length(idx),1); end
C{2,1} = param.mu*[ee; ee];

%% indexes for blocks Li
[Iall,Jall] = find(triu(ones(p)));
tmp = [Iall,Jall];
m2 = size(tmp,1);
Icomp = tmp(:,1); Jcomp = tmp(:,2);
IJ= (Icomp-1)*p+Jall;

[Ialli,Jalli] = find(triu(ones(k)));
IJi= (Ialli-1)*k+Jalli;

OPTIONS.smoothing  = 1;
OPTIONS.scale_data = 0; %% or 2;
OPTIONS.plotyes    = 0;
OPTIONS.tol        = 1e-10;

blk0=blk;
C0=C;
At0=At;

c=true;
i=0;
M=speye(p);
while c
    i = i+1;
    
    %% solve problem P_S
    
    if ~isempty(ActiveSet.I) % now we already have all the Eij+Eji atoms
        m=length(ActiveSet.I);
        if size(blk,1)>2
            blk=blk0;
            C=C0;
            At=At0;
        end
        for mm=1:m
            blk{2+mm,1} = 's'; blk{2+mm,2} = k;
            phi=sparse(ActiveSet.I{mm},(1:k)',1,p,k);
            A=kron(phi,phi);
            At{2+mm,1}  = A(IJ,IJi)';
            C{2+mm,1}   = param.lambda*speye(k,k);
        end                
        
        eta = [1; zeros(m+1,1)];
        [obj,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,eta,OPTIONS);
        M=X{1};
        obj_sup = [obj_sup obj(1)];
        
        %Remove useless blocks
        for mm=1:m
            if trace(X{2+mm})<1e-3
                ActiveSet.I{mm}=[];
                ActiveSet.k{mm}=[];
                ActiveSet.I=ActiveSet.I(~cellfun('isempty',ActiveSet.I));  
                ActiveSet.k=ActiveSet.k(~cellfun('isempty',ActiveSet.I));  
            end
        end
    end
    
    %% Cleaning, proximal steps
    
    %% get a new descent direction using truncated power iteration
    
    H = inv(M)-Sigma;

    
    if param.verbose==1
        fprintf('%d/%d   \n',i,param.max_nb_main_loop);
    end
    
    [u, kBest,val] = lmo_spsd_TPower(-H,param);
    cf=min(param.cardfun(kBest:end));

    if val<0
        currI=[];
        fprintf('   all eigs are negative\n')
        %         keyboard;
    else
        param.k=kBest;
        currI = find(u);
    end
    
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
        if ~isempty(ActiveSet.I)
            takenI= isInCell(currI,ActiveSet.I,cell2mat(ActiveSet.k)) ;
        else
            takenI=false;
        end
    end
    
    hist.varIJ=[hist.varIJ varIJ];
    flag.var=varIJ;
    
    if param.verbose==1
        fprintf(' maxIJ = %2.4e, thresh = %2.4e\n',full(maxIJ), param.mu*(1+param.epsStop));
        fprintf(' varIJ = %2.4e, thresh = %2.4e\n',full(varIJ), param.lambda*(1+param.epsStop / kBest)* param.cardfun(kBest));
        fprintf(' length(currI)=%d\n', length(currI));
    end
    
    
    if varIJ < param.lambda*cf*(1+param.epsStop) && maxIJ < param.mu*(1+param.epsStop)
        c=0;
        c=1;
%     elseif ActiveSet.atom_count>=param.max_nb_atoms
%         c=0;
    elseif takenI
        fprintf(' This support has already been added. Stopping\n');
        %c=0;
    elseif varIJ > param.lambda*cf*(1+param.epsStop)
        ActiveSet.I = [ActiveSet.I, currI];
        ActiveSet.k = [ActiveSet.k , kBest];
    else
        %c = 0;
    end
    c = i<param.max_nb_main_loop & c;
end


if param.debug==1
    if i>=param.max_nb_main_loop
        fprintf('\n max number of main loop iterations reached\n');
    end
end

%%
L={};
U={};
S=M;
m=length(ActiveSet.I);
for mm=1:m
    I= ActiveSet.I{mm};
    L{mm}=zeros(p);
    L{mm}(I,I)=X{2+mm};
    [V,D] = eig(X{2+mm});
    d=diag(D);
    ll=sum(d>1e-3);
    u = bsxfun(@times,sqrt(d(d>1e-3))',V(:,d>1e-3));
    U{mm}=zeros(p,ll);
    U{mm}(I,:)=u;
    S=S+L{mm};
end
S(abs(S)<1e-3)=0;
hist.obj_sup = obj_sup;


end





