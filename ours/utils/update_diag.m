function D = update_diag(param,inputData,L,S,D)

Z=S+D+L;
p=size(Z,1);

switch param.f
    case 1 % prox
        H = Z - inputData.Y;
        L = 1;
    case 4 % bilinear
        H = inputData.X1'*(inputData.X1*Z*inputData.X2 - inputData.Y)*inputData.X2';
        L=norm(inputData.X1,'fro')^4;
    case 5 % score matching 
%         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
        L=norm(inputData.X1,'fro');
end

%projected gradient
D=D-1/L*H;
D=spdiags(D,0,p,p);
D(D>0)=0;


end