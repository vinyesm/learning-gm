function [ Sp ] = proj_psd( S )
    S = .5*(S+S');
    [u,s] = eig(S);
    Sp = u*(s.*(s>0))*u';
end

