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

paramRL1.eps=1e-3;
paramRL1.mu=param.lambda;

paramRO.k=param.k;
paramRO.eps=1e-3;
paramRO.lambda=param.lambda;
paramRO.maxIter=1000;
paramRO.maxNbAtoms=100;

iter = 0;

% memory allocation
grad1 = zeros(p,p);
grad2 = zeros(p,p);
Y=inputData.Y;
Y2=inputData.Y2;

% main loop
tic
while iter<=param.maxIter
    
    iter = iter+1;
    
    % S sparse update
    inputData.Y=-inputData.X*init.M*inputData.X-Y;
    inputData.Y2=inputData.X*inputData.Y*inputData.X;
    [init.S, nb_iter, dg] = regL1(paramRL1,inputData,init,1);
    % hist S and update of inputData for L update
    inputData.Y=-inputData.X*init.S*inputData.X-Y;
    inputData.Y2=inputData.X*inputData.Y*inputData.X;
    [grad1, grad2]=varFenchel(inputData, init, grad1, grad2, 1);
    [~, valL]= dualOmega(-grad2,inf,param.k);
    valS = max(abs(-grad2(:)));
    dg = dualityGapSL(init, grad1, grad2, valL, valS, param);
    hist.dualityGap=[hist.dualityGap dg];
    hist.loss= [hist.loss loss(grad1,1)];
    hist.omega= [hist.omega omega(init,1)];
    hist.l1= [hist.l1 sum(abs(init.S(:)))];
    hist.objective=[hist.objective hist.loss(end)+param.lambda*hist.omega(end)+param.mu*hist.l1(end)];
    hist.time=[hist.time toc];
    hist.reldg= hist.dualityGap./hist.objective;
    
    if hist.reldg(end)<param.epsStop
        break
    end
    
    % L update
    nba = size(init.Xatoms_u,2);
    for ia=1:nba
        xatom=init.Xatoms_u(:,ia);
        qp.b = -xatom'*inputData.Y*xatom + param.lambda;
    end
    [ init, histL, qp ] = regOmegaRestricted( inputData, paramRO, set, init, qp );
    
    % hist L
    nbit=length(histL.dualityGap);
    hist.dualityGap=[hist.dualityGap histL.dualityGap];
    hist.loss= [hist.loss histL.loss];
    hist.omega= [hist.omega histL.omega];
    hist.l1= [hist.l1 sum(abs(init.S(:)))*ones(1,nbit)];
    hist.objective=[hist.objective histL.objective+param.mu*hist.l1];
    hist.time=[hist.time histL.time];
    hist.reldg= hist.dualityGap./hist.objective;
    
    if hist.reldg(end)<param.epsStop
        break
    end
    
%     keyboard;
    
end

init.M=0;
for ia=1:init.atomCount
    atom = init.atoms_u(:,ia);
    init.M = init.M + init.coeff(ia)*(atom*atom');
end

if iter>=param.maxIter
    warning('Maximum number of iterations reached. Breaking..')
end

output=init;



