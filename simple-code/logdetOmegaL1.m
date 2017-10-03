function [S,M,L,U,hist,setout] = logdetOmegaL1(Sigma,param,set)

% This code uses logdetPPA solver implemented by Professor Kim-Chuan Toh
% Chengjing Wang, Defeng Sun and Kim-Chuan Toh,
% "Solving log-determinant optimization problems by a Newton-CG primal proximal point algorithm",
% SIAM J. Optimization, 20 (2010), pp. 2994--3013
% http://www.math.nus.edu.sg/~mattohkc/LogdetPPA.html

%
% Solves problem :  min_M -likelihood(S-M; Sigma)+ lambda*Omega_S(M) + mu*|S|_1 (P)
% with likelihood(K; Sigma) := logdet(K)-tr(K*Sigma) and
% Omega_B is the norm Omega restricted to the supports contained in
% set S. M is by construction M=sum_i coeff_i*u_i*u_i'

% INPUT :
%   - Sigma : a PSD matrix X of size p*p
%   - param.lambda : lambda, regularization parameter
%   - param.k : sparsity parameter of Omega
%   - param.epsObj : stopping criterion on objective
%   - param.maxIter : maximum number of iterations
%   - param.maxNbBlocks : maximum number of atoms
%   - set : a boolean sparse matrix of size p*number_of_supports_in_S where
%   each column corresponds to a support (1 if index is in the support, 0
%   otherwise)
%   - init.M : initial solution M=sum_i coeff_i*u_i*u_i'
%   - init.S : initial solution of S, symmetric sparse matrix of dimension p*p
%   - init.coeff : coeffitients
%   - init.atoms_u : a sparse matrix of size p*param.maxNbAtoms where each
%   column is an atom u_i
%   - init.atomCount : number of active atoms
%   - init.Xatoms_u : a sparse matrix of size p*param.maxNbAtoms contains
%   X*u_i (for computational purposes)
%   - qp : when the atoms u_i are fixed the optimization on coeffs is a
%   quadratic problem quadratic problem min_coeff .5*coeff'*H*coeff+b'*coeff
%   qp.H isthe Hessian and qp.b the linear term of the quadratic problem
%
%   OUTPUT :
%   output is a structure containing the optimal M and atoms
%   hist : contains history of objective (P)
%
%   TODO : extend to multiple levels of sparsity
%
% Marina Vinyes - Ecole des Ponts ParisTech, 2017

%% init
hist.objective = [];
hist.time = [];

if set == inf
    set = [];
end
param.epsStop=1e-8;
param.debug=0;


%% set up SDP data in SDPT3 format
p=size(Sigma,1);
k=param.k;

% blocks sparse and trace
p2 = p*(p+1)/2;
b = zeros(p2,1);
C{1} = Sigma;
blk{1,1} = 's'; blk{1,2} = p;
[Iall,Jall] = find(triu(ones(p)));
tmp = [Iall,Jall];
m2 = size(tmp,1);
Icomp = tmp(:,1); Jcomp = tmp(:,2);
Itmp  = Icomp + Jcomp.*(Jcomp-1)/2;
Atmp  = spconvert([Itmp,(1:m2)',ones(m2,1); p2,m2,0]);
At{1} = Atmp;

blk{2,1} = 'l'; blk{2,2} = 2*p2;
Identity = speye(p2);
At{2,1} = [-Identity,Identity]';
idx = find(Icomp == Jcomp);
ee  = sqrt(2)*ones(m2,1);
if ~isempty(idx);
    ee(idx) = ones(length(idx),1);
end
C{2,1} = param.mu*[ee; ee];

% indexes for blocks Li
[Iall,Jall] = find(triu(ones(p)));
tmp = [Iall,Jall];
m2 = size(tmp,1);
Icomp = tmp(:,1); Jcomp = tmp(:,2);
IJ= (Icomp-1)*p+Jall;

[Ialli,Jalli] = find(triu(ones(k)));
IJi= (Ialli-1)*k+Jalli;

% param
OPTIONS.smoothing  = 1;
OPTIONS.scale_data = 0; %% or 2;
OPTIONS.plotyes    = 0;
OPTIONS.tol        = 1e-10; %% 1e-12
OPTIONS.printlevel = 0;

blk0=blk;
C0=C;
At0=At;

%% init
calls=0;
iter=0;
relobj=inf;
X{1}=speye(p);

tic
%qs= 10.^(4:-1:1);
for q=1;    
    while iter<=param.maxIter 
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------');
    fprintf('\n Call %d\n', calls);
    fprintf('-------------------------------------------------------------');
    fprintf('-------------------------------------------------------------');
        iter = iter+1;
        nblocks=size(set,2);
        if size(blk,1)>2
            blk=blk0;
            C=C0;
            At=At0;
        end
        %         for mm=1:m
        nset = size(set,2);
        for mm=1:nset
            I=find(set(:,mm));
            blk{2+mm,1} = 's'; blk{2+mm,2} = k;
            phi=sparse(I,(1:k)',1,p,k);
            A=kron(phi,phi);
            At{2+mm,1}  = A(IJ,IJi)';
            C{2+mm,1}   = param.lambda*speye(k,k);
        end
        
        eta = [1; zeros(nblocks+1,1)];
        [obj,X,y,Z,info,runhist] = logdetPPA(blk,At,C,b,eta,OPTIONS);
        calls = calls+1;
        hist.objective = [hist.objective obj(1)];
        hist.time=[hist.time toc];
        
        if length(hist.objective)>1
            relobj = (hist.objective(calls-1)-hist.objective(calls))/abs(hist.objective(calls));
        end
        
        if  relobj < param.epsObj
            fprintf('Objective stagnating. Stopping.\n')
            break;
        end
        
        %Remove useless blocks
        ii=0;
        nset=0;
        for mm=nblocks:-1:1
            ii = ii+1;
            if trace(X{2+mm})<1e-6
                set(:,mm)=[];
%                 keyboard
            else
                nset = nset+1;
                X{2+ii} = X{2+mm};
            end
        end
        
        
        %% get a new descent direction using truncated power iteration
        
        grad = inv(X{1})-Sigma;
        [newAtom, val]= dualOmega(grad,inf,param.k);
        
        %adding new block
        if val > param.lambda
            if ~isempty(set) && any(prod(bsxfun(@times, newAtom~=0, set))==1)
                warning('Repeated block found. Breaking..')
%                 break
            else
                set = [set, newAtom~=0];
            end
        else
%             keyboard;
            warning('No new block found. Breaking..')
%             break;
        end
        
    end
end


if param.debug==1
    if iter>=param.maxIter
        fprintf('\n max number of main loop iterations reached\n');
    end
end

%%
L={};
U={};
S=X{1};
M=zeros(p);
%m=size(set,2);
thresh1 = 1e-3;
thresh2 = 1e-6;
it = 0;
setout = [];
for mm=1:nset
    I=find(set(:,mm));
    [V,D] = eig(X{2+mm});
    d=diag(D);
    ll=sum(d>thresh1);
    if ll > 0
        setout = [setout set(:,mm)];
        it=it+1;
        L{it}=zeros(p);
        L{it}(I,I)=X{2+it};
        u = bsxfun(@times,sqrt(d(d>thresh1))',V(:,d>thresh1));
        U{it}=zeros(p,ll);
        U{it}(I,:)=u;
        S=S+L{it};
        M=M+L{it};
    end
end
S(abs(S)<thresh2)=0;

hist.nbcalls = calls;

end





