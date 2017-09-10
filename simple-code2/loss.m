function [val] = loss(grad1, fid)

% fid is the loss id
%           f(M)=0.5*|X*M*X-Y|_F^2 (1)
%

if fid == 1    
    val = 0.5*norm(grad1,'fro')^2;    
else
    error('this loss is not implemented');
end