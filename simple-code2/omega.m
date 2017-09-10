function [val] = omega(init, fid)

% fid is the loss id
% Omega(M) = min_c { sum_i ci | M = sum ci ui*ui'
%                   s.t ci>0 and |u_i|_0<=k }
%

if fid == 1    
    val = sum(init.coeff(1:init.atomCount));   
else
    error('this norm is not implemented');
end