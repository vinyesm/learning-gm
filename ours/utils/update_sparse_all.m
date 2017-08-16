function [S] = update_sparse_all(param,inputData,L,S)
debug=1;
spams=1;
Sold=S;
Z=S+L;
p=size(Z,1);
Y=inputData.Y-inputData.X1*L*inputData.X2;
X=inputData.X1;


if spams
    param_spams.loss='cur';
    param_spams.lambda=param.mu;
    param_spams.regul='l1';
    param_spams.max_it=1000;
    param_spams.verbose=1;
    param_spams.tol=param.epsStop;
%     .5*norm(Y-X*Sold*X,'fro')^2+param.mu*sum(abs(Sold(:)));
    [W optim]=mexFistaFlat(Y,X,Sold,param_spams);
%     .5*norm(Y-X*W*X,'fro')^2+param.mu*sum(abs(W(:)));
    S=W;
    % keyboard;
end

if ~spams
switch param.f
    case 1 % prox
        %         H = Z - inputData.Y;
        Lip = 1;
    case 4 % bilinear
        loss= .5*norm((inputData.X1*Z*inputData.X2 - inputData.Y),'fro')^2;
        obj0= loss + param.mu*sum(abs(S(:)));
        grad_S = inputData.X1'*(inputData.X1*Z*inputData.X2 - inputData.Y)*inputData.X2';
        Lip=norm(inputData.X1,'fro')^4;
    case 5 % score matching
        %         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        %         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        Lip=norm(inputData.X1,'fro');
end

% Proximal gradient method
t=1;
Snew=S-t*grad_S;
Snew=wthresh(Snew,'s',t*param.mu);

% keyboard;
while loss- 0.5*norm((inputData.X1*(Snew+L)*inputData.X2 - inputData.Y),'fro')^2 + trace(grad_S'*(Snew-S))+1/(2*t)*norm(Snew-S,'fro')^2<0
    t=t/2;
    Snew=S-t*grad_S;
    Snew=wthresh(Snew,'s',t*param.mu);
end
S=Snew;

if debug
    obj1=.5*norm((inputData.X1*(Snew+L)*inputData.X2 - inputData.Y),'fro')^2+param.mu*sum(abs(Snew(:)));
    if obj0<=obj1
        fprintf('error update S : objective does not decrease\n');
        keyboard;
    end
end

%computing dg
% M=L+S;
% H = inputData.X1'*(inputData.X1*(S+L)*inputData.X2 - inputData.Y)*inputData.X2';
% temp_l1 = max(abs(H(:)));
% shrink=min(1,param.mu / temp_l1);
% kappa = shrink * (inputData.X1*M*inputData.X2-inputData.Y);
% G=inputData.X1*kappa*inputData.X2;
% dg_f=.5*norm(inputData.X1*M*inputData.X2 - kappa -inputData.Y,'fro')^2;
% dg_S=param.mu*sum(abs(S(:)))+trace(G*S);
% dg=dg_f+dg_S;
% 
% if dg_f<0 || dg_S<0  abs(dg)>1e-10
%     fprintf('Negative duality gap\n');
%     keyboard;
%     %     error('Negative duality gap=%f, gapLoss=%f gapPen=%f\n',dualityGap, gapLoss, gapPen);
%     dg=abs(dg);
% end

end


