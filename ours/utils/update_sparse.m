function S = update_sparse(param,inputData,L,S,D)
debug=1;
Sold=S;
Z=S+D+L;
p=size(Z,1);

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
% Snew=S;
% Sigma=inputData.X1'*inputData.X1;
% R= inputData.X1'*(inputData.X1*(D+L)*inputData.X2 - inputData.Y)*inputData.X2';
I=speye(p)==1;
% t=2/Lip;
t=1;
Snew=S-t*grad_S;
Snew(I)=0;
Snew=wthresh(Snew,'s',t*param.mu);

% keyboard;
while loss- 0.5*norm((inputData.X1*(Snew+D+L)*inputData.X2 - inputData.Y),'fro')^2 + trace(grad_S'*(Snew-S))+1/(2*t)*norm(Snew-S,'fro')^2<0
    t=t/2;
    Snew=S-t*grad_S;
    Snew(I)=0;
%     Snew=(1-eye(p)).*Snew;
%     keyboard;
    Snew=wthresh(Snew,'s',t*param.mu);
%     obj1=.5*norm((inputData.X1*(Snew+L+D)*inputData.X2 - inputData.Y),'fro')^2 +param.mu*sum(abs(Snew(:)));
%     grad_S=Sigma*Snew*Sigma+R;
%
%     t=t/2;
end
% keyboard;
S=Snew;
% S=(Snew-param.mu*(ones(p)-eye(p))).*sign(Snew);
% S= wthresh(Snew,'s',param.mu);

% keyboard;
% 
% S((speye(p)>0))=0;

if debug
    obj1=.5*norm((inputData.X1*(Snew+L+D)*inputData.X2 - inputData.Y),'fro')^2+param.mu*sum(abs(Snew(:)));
    if obj0<=obj1
        fprintf('error update S : objective does not decrease\n');
        keyboard;
    end
%     keyboard; % [obj0 obj1]
end