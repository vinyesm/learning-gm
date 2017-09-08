function [S, nb_iter, reldg] = regL1(param,inputData,init,fid)

%relative duality gap dg/objective

if fid==1
    param_spams.loss='cur';
    param_spams.lambda=param.mu;
    param_spams.regul='l1';
    param_spams.max_it=10000;
    param_spams.verbose=1;
    param_spams.tol=param.eps;
    [S optim]=mexFistaFlat(inputData.Y,inputData.X,full(init.S),param_spams);
    S=sparse(S);
    reldg=optim(3);
    %dg=optim(3)*optim(1);
    nb_iter=optim(4);
    if nb_iter>=param_spams.max_it
    end
end



end


