function [obj, loss, pen] = get_val_l1_omega_02(L,S,D,inputData,param,ActiveSet)

Z=S+D+L;

cf=zeros(ActiveSet.atom_count,min(ActiveSet.atom_count,1));
card=sum(ActiveSet.atoms(:,1:ActiveSet.atom_count)~=0);
for i=1:ActiveSet.atom_count
    cf(i)=min(param.cardfun(card(i):end));
end

pen = param.lambda*dot(cf,ActiveSet.alpha)+param.mu*(.5*sum(abs(S(:))));

switch param.f
    case 1 % prox
        loss = norm(Z - inputData.Y,'fro')^2;
    case 4 % bilinear
        loss = norm(inputData.X1*Z*inputData.X2 - inputData.Y,'fro')^2;
    case 5 % score matching 
        loss = .5*trace(inputData.X*Z*Z)+trace(Z);
end

obj=pen+loss;

end