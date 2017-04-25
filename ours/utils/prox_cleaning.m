function  [Z2,ActiveSet]=prox_cleaning(Z1,Z2,S,ActiveSet,param)

S05=S^.5;
Z=Z2;
p=size(ActiveSet.atoms,1);
nb=length(ActiveSet.I);
T=10;
L= norm(S,'fro');%lipshitz

J=false(1,ActiveSet.atom_count);
%building blocks
for i=1:nb
    supp=ActiveSet.I{i};
    len=length(supp);
    block=zeros(len);
    K=false(1,ActiveSet.atom_count);
    for at=1:ActiveSet.atom_count
        u=ActiveSet.atoms(:,at);
        %         keyboard;
        if len==sum(u~=0) && all(supp==find(u~=0))
            K(at)=1;
            block=block+ActiveSet.alpha(at)*(u*u');
        end
    end
end
%proximal steps
for t=1:T
    order=randperm(nb); 
    for i=order
        block=ActiveSet.block{i};
        Z2i=Z-block;
        grad=S*(Z1+Z)*S+S;
        supp=sum(block(1,:)~=0);
        cf=param.cardfun(supp);
%         ActiveSet.block{i}=prox_trace(block-grad/L,param.lambda*cf/L);
    end    
end


end
