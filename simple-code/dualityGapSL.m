function dg = dualityGapSL(init, grad1, grad2, valL, valS, param)

if valL<0
    shrink=min(1,param.mu / valS);
else
    shrink=min(min(1, param.lambda / valL),param.mu / valS);
end

% duality gap
dgf=.5*(1-shrink)*norm(grad1,'fro')^2;
dgS=param.mu*sum(abs(init.S(:)))+shrink*trace(init.S*grad2);
dgL= param.lambda*sum(init.coeff(1:init.atomCount)) + shrink*trace(init.M*grad2);
dg = dgf+dgL+dgS;

if param.verbose>=2
    fprintf('dgf=%f  dgL=%f dgS=%f dg_global=%f  (<1)\n',dgf,dgL,dgS,dg);
end

