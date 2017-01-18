function [ M ] = rand_sym_mat( sz,rank )

spec=[(rand(1,rank)-.5) zeros(1,sz-rank)];
[Q,R]=qr(randn(sz));
M=Q'*diag(spec)*Q;
M=(M+M')/2;

end

