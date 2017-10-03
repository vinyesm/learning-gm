function [order]=grayorder(M)
[n,p]=size(M);

x=2.^(0:(p-1))';
b=M*x;

c=zeros(n,1);


for i=1:n
    c(i)=gray2int(b(i));
end

[res order]=sort(c);


end