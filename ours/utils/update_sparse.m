function S = update_sparse(param,L,S,D)

switch param.f
%     case 1 % prox
%         H = Z - inputData.Y;

%     case 4 % bilinear
%         H = inputData.X1'*(inputData.X1*Z*inputData.X2 - inputData.Y)*inputData.X2';
%     case 5 % score matching 
% %         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
%         H = .5*(inputData.X*Z+Z*inputData.X)-inputData.Y;
end


end