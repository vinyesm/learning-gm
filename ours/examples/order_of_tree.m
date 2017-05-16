function [C,I]=order_of_tree(Z)
n=size(Z,1)+1;
C=NaN(2*n-1,2*n-1);
for i=(n-1):-1:1
    k=sum(~isnan(C(n+i,:)));
    s1=Z(i,1);
    s2=Z(i,2);
    C(s1,1:k)=C(n+i,1:k);
    C(s2,1:k)=C(n+i,1:k);
    C(s1,k+1)=-1;
    C(s2,k+1)=1;
%     C2=C;
%     C2(isnan(C2))=0;
%     imagesc(C2);
%     pause();
end

k=sum(~all(isnan(C),1));
C=C(1:n,1:k);
C(isnan(C))=0;
x=2.^(k-1:-1:0)';
v=C*x;
[~,I]=sort(v);