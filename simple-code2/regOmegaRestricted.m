function [ output, hist, qp] = regOmegaRestricted( inputData, param, set, init, qp )
% Solves problem :  min_M |X*M*X-Y|_F^2 + lambda*Omega_S(M)  (P)
% where Omega_B is the norm Omega restricted to the supports contained in
% set S. Where M is bay construction M=sum_i coeff_i*u_i*u_i'
%
% INPUT :
%   - inputData.X : a PSD matrix X of size p*p
%   - inputData.Y : Y of size p*p
%   - param.lambda : lambda, regularization parameter
%   - param.k : sparsity parameter of Omega
%   - param.eps : stopping criterion on the relative duality gap
%   - param.maxIter : maximum number of iterations
%   - param.maxNbAtoms : maximum number of atoms
%   - set : a boolean sparse matrix of size p*number_of_supports_in_S where
%   each column corresponds to a support (1 if index is in the support, 0
%   otherwise)
%   - init.M : initial solution M=sum_i coeff_i*u_i*u_i'
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
%   hist : contains history of objective, loss, penalty and dualityGap of(P)
%
%   TODO : extend to multiple levels of sparsity
%
% Marina Vinyes - Ecole des Ponts ParisTech, 2017

if nargin < 3
    error('Not enough input arguments');
end

p = size(inputData.Y,1);
paramQP.pbId=1;

if ~isfield(inputData,'X2')
    inputData.X2=inputData.X*inputData.X;
    inputData.Y2=inputData.X*inputData.Y*inputData.X;
end

if ~isfield(param,'verbose')
    param.verbose=1;
end

if nargin < 4
    init.coeff = sparse(param.maxNbAtoms,1);
    init.atoms_u = sparse(p,param.maxNbAtoms);
    init.Xatoms_u = sparse(p,param.maxNbAtoms);
    init.X2atoms_u = sparse(p,param.maxNbAtoms);
    init.atomCount=0;
    %build qp.H and qp.b
    [ qp ] = buildQP( inputData, init.Xatoms_u(:,1:init.atomCount), param.lambda, paramQP );
elseif nargin < 5
    % check consistency init and qp
end

hist.objective = zeros(1,param.maxIter);
hist.time = zeros(1,param.maxIter);
hist.dualityGap = zeros(1,param.maxIter);
hist.loss = zeros(1,param.maxIter);
hist.omega = zeros(1,param.maxIter);
hist.reldgl1 = zeros(1,param.maxIter);
hist.reldg = zeros(1,param.maxIter); %GO

paramAS.max_iter=1e3;
paramAS.epsilon=1e-14;
paramAS.debug_mode=false;
paramAS.ws=true;

count = 0;
sloppyCount = 0;

% memory allocation
grad1 = zeros(p,p);
grad2 = zeros(p,p);
if isfield(init,'S')
    l1 = sum(abs(init.S(:)));
end

% main loop
firstPass=false; %%GO: I changed this to false, I do not understand its usefulness
new_atom_added = false;
tic
while count<=param.maxIter
    count = count+1;
    if init.atomCount>param.maxNbAtoms
        warning('Maximum number of atoms reached. Breaking..')
        keyboard
        break;
    end
    
    if ~firstPass
        %active-set
        if init.atomCount>0
            c0 = [init.coeff(1:init.atomCount-1); 0];
            %         if abs(c0'*qp.H*c0-norm(inputData.X*init.M*inputData.X, 'fro')^2)/abs(c0'*qp.H*c0)>1e-14
            %             keyboard;
            %         end
            %         if abs(-trace(inputData.X*init.M*inputData.X*inputData.Y)+param.lambda*sum(c0)-dot(c0,qp.b))>1e-14
            %             keyboard;
            %         end
            %         g = qp.H*c0+qp.b;
            %         if g(end)>0
            %             keyboard
            %         end
            if ~new_atom_added,
                idx_atom=-1;
            end
            [c,A,~,ng]=asqp2(qp.H,-qp.b,c0,paramAS,new_atom_added,idx_atom,init.atoms_u(:,init.atomCount));
            K = c>0;
            qp.H=qp.H(K,K);
            qp.b=qp.b(K);
            %fprintf('TODO : if Y changes qp.b has to be updated at each iteration\n');
            init.atoms_u(:,1:sum(K)) = init.atoms_u(:,K);
            init.Xatoms_u(:,1:sum(K)) = init.Xatoms_u(:,K);
            init.X2atoms_u(:,1:sum(K)) = init.X2atoms_u(:,K);
            init.atomCount = sum(K);
            init.coeff(1:sum(K))= c(K);
            init.coeff(sum(K)+1:end)= 0;
            init.atoms_u(:,sum(K)+1:end) = 0;
            init.Xatoms_u(:,sum(K)+1:end) = 0;
            init.X2atoms_u(:,sum(K)+1:end) = 0;
        end
    end
    
    %new atom
    [grad1, grad2] = varFenchel(inputData, init, grad1, grad2, 1);
    [newAtom, val]= dualOmega(-grad2,set,param.k);
    %     [newAtom, val]= dualOmega(-grad2,inf,param.k);
    
    %test adding same atom
    if init.atomCount>0
        scal = sum(bsxfun(@times,newAtom,init.atoms_u(:,1:init.atomCount)));
        [scalmin, idx]=max(abs(scal));
        if scalmin>.9
            %fprintf('adding same atom\n');
            %             keyboard
        end
    end
    
    %dualityGap
    if mod(count,10)==0
        sloppyCount = sloppyCount+1;
        init.M = 0;
        for ia=1:init.atomCount
            atom = init.atoms_u(:,ia);
            init.M = init.M + init.coeff(ia)*(atom*atom');
        end
        dg = dualityGapL(init, grad1, grad2, val, param);
        hist.dualityGap(sloppyCount)=dg;
        hist.loss(sloppyCount)=loss(grad1,1);
        hist.omega(sloppyCount)=omega(init,1);
        hist.objective(sloppyCount)=hist.loss(sloppyCount)+param.lambda*hist.omega(sloppyCount);
        hist.reldg(sloppyCount)=hist.dualityGap(sloppyCount)/hist.objective(sloppyCount);
        if isfield(init,'S')
            valS = max(abs(-grad2(:)));
            dgl1 = dualityGapS(init, grad1, grad2, valS, param);
            hist.reldgl1(sloppyCount)=dgl1/(hist.loss(sloppyCount)+param.mu*l1); 
        end
        hist.time(sloppyCount)=toc;
        if dg/hist.objective(sloppyCount)<param.eps
            break;
        end
    end
    
    %adding new atom
    if val > param.lambda
        newXatom=inputData.X*newAtom;
        newX2atom=inputData.X2*newAtom;
        init.atoms_u(:,init.atomCount+1) = newAtom;
        init.Xatoms_u(:,init.atomCount+1) = newXatom;
        init.X2atoms_u(:,init.atomCount+1) = newX2atom;
        [ qp ] = updateQP( inputData, qp, init.Xatoms_u(:,1:init.atomCount), newXatom, param.lambda, paramQP );
        init.atomCount = init.atomCount+1;
        new_atom_added = true;
        idx_atom = init.atomCount;
    else
        warning('No new atom found. Breaking..')
             %keyboard;
        break;
    end
    
    firstPass=false;
end

if count>=param.maxIter
    warning('Maximum number of iterations reached. Breaking..')
end

hist.dualityGap = hist.dualityGap(1:sloppyCount);
hist.loss = hist.loss(1:sloppyCount);
hist.omega = hist.omega(1:sloppyCount);
hist.objective = hist.objective(1:sloppyCount);
hist.reldgl1 = hist.reldgl1(1:sloppyCount);
hist.reldg = hist.reldg(1:sloppyCount);%GO
hist.time = hist.time(1:sloppyCount);
output=init;
