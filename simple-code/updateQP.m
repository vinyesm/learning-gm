function [ qp ] = updateQP( inputData, qp, Xatoms_u, newXatom, lambda, param)
% Builds the corresponding the Hessian and the linear term (Q,b)
%   min_x 0.5*x*H*x + b'*x (QP)
% by updating the old Hessian where X hasn't changed,
% for different problems
% pbId :
% min_x |X*M*X-Y|_F^2 + lambda*sum_i x_i with M=sum_i x_i*u_i*u_i'   (1)
%
% INPUTS :
% - param.pbId is an integer designing the considered problem
% - inputData contains X and Y
% - Xatoms_u : a sparse matrix of size p*param.maxNbAtoms contains
%   X*u_i (for computational purposes)
% - newXatom : the added atom, multiplied by X

nba = size(Xatoms_u,2);

if param.pbId==1
    qpNew.H=zeros(nba+1,nba+1);
    if nba>0
        qpNew.H(1:nba,1:nba)=qp.H;
        v=Xatoms_u(:,1:nba)'*newXatom;
        v=v.^2;
        qpNew.H(nba+1,1:nba)=v';
        qpNew.H(1:nba,nba+1)=v;
    end
    qpNew.H(nba+1,nba+1)=dot(newXatom,newXatom)^2;
    qpNew.b = zeros(nba+1,1);
    for ia=1:nba
        xatom=Xatoms_u(:,ia);
        qpNew.b(ia) = -xatom'*inputData.Y*xatom + lambda;
    end
    qpNew.b(nba+1) = -newXatom'*inputData.Y*newXatom + lambda;
    qp=qpNew;
else
    error('This problem Id is not considered');
end

end

