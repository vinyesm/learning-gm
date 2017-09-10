function dg = dualityGapL(init, grad1, grad2, val, param)

if val<0
    shrink=1;
else
    shrink=min(1, param.lambda / val);
end

% duality gap
dgf=.5*(1-shrink)*norm(grad1,'fro')^2;
dgL= param.lambda*sum(init.coeff(1:init.atomCount)) + shrink*trace(init.M*grad2);
dg = dgf+dgL;
if param.verbose>=2
    fprintf('dgf=%f  dgL=%f  dg_all_L=%f  (<1)\n',dgf,dgL,dg);
end