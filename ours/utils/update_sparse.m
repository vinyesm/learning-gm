function S = update_sparse(param,inputData,L,S,D)
debug=1;
Z=S+D+L;
p=size(Z,1);

switch param.f
    case 1 % prox
        %         H = Z - inputData.Y;
        Lip = 1;
    case 4 % bilinear
        obj0=.5*norm((inputData.X1*Z*inputData.X2 - inputData.Y),'fro')^2;
        grad_S = inputData.X1'*(inputData.X1*Z*inputData.X2 - inputData.Y)*inputData.X2';
        Lip=norm(inputData.X1,'fro')^2;
    case 5 % score matching
        %         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        %         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        Lip=norm(inputData.X1,'fro');
end

%prox l1
I=speye(p)==0;
S(I)=S(I)-1/Lip*grad_S(I);
S= wthresh(S,'s',param.mu);

% S((speye(p)>0))=0;

if debug
    obj1=.5*norm((inputData.X1*(S+L+D)*inputData.X2 - inputData.Y),'fro')^2;
    if obj0<=obj1
        fprintf('error : objective does not decrease\n');
        keyboard;
    end
%     keyboard; % [obj0 obj1]
end