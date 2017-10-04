function dg = dualityGapS(init, grad1, grad2, val, param)

shrink=min(1, param.mu / val);


% duality gap
dgf=.5*(1-shrink)*norm(grad1,'fro')^2;
dgS= param.mu*sum(abs(init.S(:))) + shrink*trace(init.S*grad2);
dg = dgf+dgS;

if dgS<0
    keyboard;
end

if param.verbose>=2
    fprintf('dgf=%f  dgS=%f  dg_all_S=%f  (<1)\n',dgf,dgS,dg);
end
