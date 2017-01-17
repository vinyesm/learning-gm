function [ Z ] = prox_omega(X0,lambda)

Z=X0-proj_omega(X0,lambda);


end

