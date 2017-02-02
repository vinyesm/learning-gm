
clc;clear all:

p=3;
S=randn(p);
S=S*S';

E=speye(p*p);
pairs=fullfact([p p]);
ii=pairs(:,1)>=pairs(:,2);
pairs=[pairs(ii,1) pairs(ii,2)];

J=1:(2*length(pairs)-p);
I=zeros(1,2*length(pairs)-p);
val=zeros(1,2*length(pairs)-p);

count=1;
for j=1:length(pairs);
    if pairs(j,1)==pairs(j,2)
        J(count)=j;
        I(count)= (pairs(j,2)-1)*p+pairs(j,1);
        val(count)=1;
        count=count+1;
    else
        J(count)=j;
        J(count+1)=j;
        I(count)= (pairs(j,2)-1)*p+pairs(j,1);
        val(count)=.5;
        I(count+1)= (pairs(j,1)-1)*p+pairs(j,2);
        val(count+1)=.5;
        count=count+2;
    end
end

atoms=sparse(I,J,val);
atoms_l1_sym=[atoms -atoms];

