function D = update_diag(param,inputData,L,S,D)

debug=1;

param_as.max_iter=1e3;
param_as.epsilon=1e-14;
param_as.debug_mode=false;
param_as.ws=true;


Z=S+D+L;
p=size(Z,1);
d0=diag(D);

switch param.f
    case 1 % prox
%         H = Z - inputData.Y;
        Lip = 1;
    case 4 % bilinear
        obj0=.5*norm((inputData.X1*Z*inputData.X2 - inputData.Y),'fro')^2;
%         H = inputData.X1'*(inputData.X1*Z*inputData.X2 - inputData.Y)*inputData.X2';
        Sigma=inputData.X1*inputData.X1;
        C2=Sigma.^2;
        grad_D= C2*diag(D)+diag(Sigma+Sigma*(S+L)*Sigma);
        Lip=norm(C2,'fro');
    case 5 % score matching 
%         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        Lip=norm(inputData.X1,'fro');
end



%%projected gradient
% d0=d0-1/Lip*grad_D;
% d=min(d0,0);
% D=spdiags(d,0,p,p);
% 
% if debug 
%     obj1=.5*norm((inputData.X1*(S+L+D)*inputData.X2 - inputData.Y),'fro')^2;
%     if obj0<obj1
%         fprintf('error update D: objective increases\n');
%         keyboard;
%     end
% end

% exact minimisation
% asqp2(Q,b,c0,param,new_atom_added,idx_atom)
% keyboard;
Hall=Sigma.^2;%ok
fall= -diag(inputData.X1*((+inputData.X1*(S+L)*inputData.X2)-inputData.Y)*inputData.X1);
c0=-d0;
[c,A,nbpivot,ng]= asqp2(Hall,-fall,c0,param_as,0,0);
D=-diag(c);

if debug 
    obj1=.5*norm((inputData.X1*(S+L+D)*inputData.X2 - inputData.Y),'fro')^2;
    if obj0<obj1
        fprintf('error update D: objective increases\n');
        keyboard;
    end
end

% keyboard;


% keyboard;
end