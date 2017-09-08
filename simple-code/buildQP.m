function [ qp ] = buildQP( inputData, Xatoms_u, lambda, param )
% Builds the corresponding the Hessian and the linear term (Q,b) 
%   min_x 0.5*x*H*x + b'*x (QP)
% for different problems 
% pbId :  
% min_x |X*M*X-Y|_F^2 + lambda*sum_i x_i with M=sum_i x_i*u_i*u_i'   (1)
% 
% INPUTS :
% - param.pbId is an integer designing the considered problem
% - inputData contains X and Y
% - Xatoms_u : a sparse matrix of size p*param.maxNbAtoms contains 
%   X*u_i (for computational purposes)

nba = size(Xatoms_u,2);

if param.pbId==1
    qp.H = Xatoms_u'*Xatoms_u;
    qp.H = (qp.H).^2;
    qp.b = zeros(nba,1);
    for ia=1:nba
        xatom=Xatoms_u(:,ia);
        qp.b = -xatom'*inputData.Y*xatom + lambda;
    end
else
    error('This problem Id is not considered');
end


end

