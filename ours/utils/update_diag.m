function D = update_diag(param,inputData,L,S,D)
debug=1;
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
        C=inputData.X1*inputData.X1;
        C2=C.^2;
        grad_D= C2*diag(D)+diag(C+C*(S+L)*C);
        Lip=norm(C2,'fro');
    case 5 % score matching 
%         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        Lip=norm(inputData.X1,'fro');
end



%%projected gradient
% d0=d0-1/Lip*grad_D;
% d=min(d0,0);
% D=spdiags(d,0,p,p);

%%projected exact minimimization
C=inputData.X1*inputData.X1;
C2=C.^2;
% d0=C2\diag(C-C*(S+L)*C); % if input.Y=eye
% d=max(d0,0);
d0=-C2\diag(C-C*(S+L)*C); % if input.Y=-eye
d=min(d0,0);
D=spdiags(d,0,p,p);

if debug 
    obj1=.5*norm((inputData.X1*(S+L+D)*inputData.X2 - inputData.Y),'fro')^2;
    if obj0<obj1
        fprintf('error update D: objective increases\n');
%         keyboard;
    end
end

% keyboard;
end