function [ output_args ] = sparse_omega_lgm( input_args )
% min_(A,M,S) .5|C^.5*A*C^.5|^2 + mu|S|_1 + lambda Omega_psd,k(M)
% s.t. A>=0 and M>=0
% using ADMM

%rho=1/L;
%init A,M,S,E

% iter
% A with ls eta=norm(G,'fro')^2/(norm(C^.5*G*C^.5,'fro')^2 + rho*norm(G,'fro')^2)
% S soft thresholding
% M = prox_omega(
% E + E+rho*(S-M-A);

end

