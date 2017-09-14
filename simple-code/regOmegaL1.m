function [ output, hist ] = regOmegaL1( inputData, param, set, init, qp )
% Solves problem :  min_M |X*(S-M)*X-Y|_F^2 + lambda*Omega_S(M) + mu*|S|_1 (P)
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
    init.S = sparse(p,p);
    init.M = zeros(p,p);
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

hist.objective = [];
hist.time = [];
hist.dualityGap = [];
hist.loss = [];
hist.omega = [];
hist.l1 = [];
hist.reldg = [];
hist.reldgl1 = [];
hist.reldgom = [];

paramRL1.mu=param.lambda;

paramRO.k=param.k;
paramRO.lambda=param.lambda;
paramRO.maxIter=1000;
paramRO.maxNbAtoms=100;
paramRO.verbose=param.verbose;
paramRO.mu=param.mu;

iter = 0;

% memory allocation
grad1 = zeros(p,p);
grad2 = zeros(p,p);
Y=inputData.Y;
Y2=inputData.Y2;

% main loop
tic

%qs= 10.^(4:-1:1);
for q=1;
    qeps=q*param.epsStop;
    paramRO.eps=10*qeps;
    paramRL1.eps=10*qeps;
    while iter<=param.maxIter
        
        iter = iter+1;
        
        % S sparse update
        inputData.Y=-inputData.X*init.M*inputData.X-Y;
        inputData.Y2=inputData.X*inputData.Y*inputData.X;
        [init.S, nb_iter, reldgl1] = regL1(paramRL1,inputData,init,1);
        
        % hist S and update of inputData for L update
%         g1 = inputData.X*(init.M+init.S)*inputData.X+Y;
%         g2 = inputData.X*(g1)*inputData.X;
        inputData.Y=-inputData.X*init.S*inputData.X-Y;
        inputData.Y2=inputData.X*inputData.Y*inputData.X;
        [grad1, grad2]=varFenchel(inputData, init, grad1, grad2, 1);
        [~, valLall]= dualOmega(-grad2,inf,param.k);
        [~, valLset]= dualOmega(-grad2,set,param.k);
        valL=max(valLall,valLset);
        valS = max(abs(-grad2(:)));
        dg = dualityGapSL(init, grad1, grad2, valL, valS, param);
        dgL= dualityGapL(init, grad1, grad2, valL, param);
        dgS= dualityGapS(init, grad1, grad2, valS, param);% to check if dgS/obj
        hist.dualityGap=[hist.dualityGap dg];
        hist.loss= [hist.loss loss(grad1,1)];
        hist.omega= [hist.omega omega(init,1)];
        hist.l1= [hist.l1 sum(abs(init.S(:)))];
        hist.objective=[hist.objective hist.loss(end)+param.lambda*hist.omega(end)+param.mu*hist.l1(end)];
        hist.reldgl1= [hist.reldgl1 reldgl1];
        hist.reldgom= [hist.reldgom dgL/(hist.loss(end)+param.lambda*hist.omega(end))];
        hist.reldg= [hist.reldg  dg./hist.objective(end)];
        hist.time=[hist.time toc];
        
        % checking if spams and ours compute the same duality gap for S
        if param.verbose>=2
            objS = hist.loss(end)+param.mu*hist.l1(end);
            fprintf('Security check : reldgl1=%f   rel_dg_allS=%f (should be equal)\n', reldgl1, dgS/objS);
            err = (reldgl1-dgS/objS)^2;
            if err>1e-5
                warning(['spams and ours return different duality gaps for S, l2_error=' num2str(err)])
            end
        end
        
        if param.verbose>=1
            fprintf('reldgom=%f  reldgl1=%f reldgglobal=%f  (<1)\n',hist.reldgom(end),hist.reldgl1(end),hist.reldg(end));
        end
        
        if hist.reldg(end)<qeps
            keyboard
            break
        end
        
        % L update 
        if param.verbose>=1
            fprintf('\n')
        end
        nba = size(init.Xatoms_u,2);
        for ia=1:nba
            xatom=init.Xatoms_u(:,ia);
            qp.b = -xatom'*inputData.Y*xatom + param.lambda;
        end
        [ init, histL, qp] = regOmegaRestricted( inputData, paramRO, set, init, qp );
        
        % hist L
        nbit=length(histL.dualityGap);
        [grad1, grad2]=varFenchel(inputData, init, grad1, grad2, 1);
        [~, valLall]= dualOmega(-grad2,inf,param.k);
        [~, valLset]= dualOmega(-grad2,set,param.k);
        valL=max(valLall,valLset);
        valS = max(abs(-grad2(:)));
        dg = dualityGapSL(init, grad1, grad2, valL, valS, param);
        hist.dualityGap=[hist.dualityGap dg];
        hist.loss= [hist.loss histL.loss];
        hist.omega= [hist.omega histL.omega];
        hist.l1= [hist.l1 sum(abs(init.S(:)))*ones(1,nbit)];
        hist.objective=[hist.objective histL.objective+param.mu*sum(abs(init.S(:)))*ones(1,nbit)];
        hist.reldgl1= [hist.reldgl1 histL.reldgl1];
        hist.reldgom= [hist.reldgom histL.reldg];
        hist.reldg= [hist.reldg  dg./hist.objective(end)];
        hist.time=[hist.time histL.time];
        
        if param.verbose>=1
            fprintf('reldgom=%f  reldgl1=%f reldgglobal=%f  (<1)\n',hist.reldgom(end),hist.reldgl1(end),hist.reldg(end));
        end
        
        if hist.reldg(end)<qeps
            keyboard
            break
        end
        
        %security break
        %         if abs(hist.obj(end-1)-hist.obj(end))<param.epsStop
        %             break
        %         end
        
    end
    
end

init.M=0;
for ia=1:init.atomCount
    atom = init.atoms_u(:,ia);
    init.M = init.M + init.coeff(ia)*(atom*atom');
end

if iter>=param.maxIter
    warning('Maximum number of iterations reached.')
end


init.atoms_u=init.atoms_u(:,1:init.atomCount);
init.atoms_Xu=init.Xatoms_u(:,1:init.atomCount);
init.atoms_X2u=init.X2atoms_u(:,1:init.atomCount);
init.coeff=init.coeff(1:init.atomCount);
output=init;



