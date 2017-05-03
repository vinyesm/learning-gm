function  [Z2,ActiveSet]=prox_cleaning(Z1,Z2,S,ActiveSet,param,T,debug)

lam=param.lambda;
debug=0;



if debug
    obj=[];
end

S05=S^.5;
Z=Z2;
p=size(Z1,1);
nb=length(ActiveSet.I);
L= norm(S,'fro')^2;%lipshitz

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
            block=block+ActiveSet.alpha(at)*(u(supp)*u(supp)');
        end
    end
    ActiveSet.block{i}=block;
    ActiveSet.supp{i}=supp;
end

%proximal steps
count=1;
for t=1:T
    order=randperm(nb); 
    %activeSet.atoms=[];
    ActiveSet.alpha=[];
    na=1;%nb atoms
    for i=order
        if debug
            reg=0;
            for j=1:nb
                cf=min(param.cardfun(size(ActiveSet.block{j},1):p));
                reg=reg+cf*trace(ActiveSet.block{j});
            end
            obj(count)=.5*norm(S05*(Z1+Z)*S05+eye(p),'fro')^2 + lam*reg;
            count=count+1;
        end
        supp=ActiveSet.supp{i};
        block=ActiveSet.block{i};
        Zi=Z(supp,supp)-block;
        grad=S*(Z1+Z)*S+S;
        grad=grad(supp,supp);
        cf=param.cardfun(size(block,1));
        % proximal step
        M=block-grad/L;
        M=.5*(M+M');
        [U,D]=eig(M);
        D=diag(D);
        ds=soft_threshold(D,lam*cf/L);
        ds(ds<0)=0;
        Us=U(:,ds>0);
        %update
        nai=sum(ds>0);
        ActiveSet.block{i}=U*diag(ds)*U';
        ActiveSet.alpha=[ActiveSet.alpha;cf*ds(ds>0)];
        ActiveSet.atoms(:,na:na+nai-1)=0;
        for ii=1:nai
            ActiveSet.atoms(:,na+ii-1)=sparse(supp,ones(length(supp),1),Us(:,ii)./sqrt(cf),p,1);
        end
        ActiveSet.atom_count=na+nai-1;
        na=na+nai;
        Z(supp,supp)=Zi+ActiveSet.block{i};
    end    
end
ActiveSet.atoms=ActiveSet.atoms(:,1:ActiveSet.atom_count);

if debug
    figure(10);clf;
    plot(obj);
    keyboard;
end


end
