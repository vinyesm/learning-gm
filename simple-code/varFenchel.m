function [grad1, grad2] = varFenchel(inputData, init, grad1, grad2, fid)
% varFenchel computes gradients used for fencel duality gap computation
% fid is the function id
%           f(M)=0.5*|X*M*X-Y|_F^2 (1)
%
debug=0;
grad1=0;
grad2=0;

if fid==1
    if debug
        grad10=inputData.X*init.M*inputData.X-inputData.Y;
        grad20=inputData.X*grad10*inputData.X;
    end
    for ia=1:init.atomCount
        Xatom = init.Xatoms_u(:,ia);
        X2atom = init.X2atoms_u(:,ia);
        grad1 = grad1 + init.coeff(ia)*(Xatom*Xatom');
        grad2 = grad2 + init.coeff(ia)*(X2atom*X2atom');
    end
    grad1=grad1-inputData.Y;
    grad2=grad2-inputData.Y2;
%     grad1=grad10;
%     grad2=grad20;
    if debug
        if norm(grad10-grad1,'fro')^2>1e-16 ||  norm(grad20-grad2,'fro')^2>1e-16 
            error('error in gradient computation\n')
        end
    end
else
    error('varFenchel not implemented for this function');
end


end