function [uBest, valBest, supportBest,val_list]= dualOmega(A,set,k)
% Computes the dual of Omega at M, where
% Omega(M) = min_c { sum_i cf(ki)*ci | M = sum ci ui*ui'
%                   s.t ci>0 and |u_i|_0=ki }
%
% If set~= inf it computes the dualof Omega restricted to the sets in set
% When set=inf an heuristic is called to compute the dual
both=0;
B=0.5*(A+A');
% emin=eigs(B,1,'sa');
emin=-norm(B,'fro')/2;

if emin<0
    B=B-2*emin*eye(size(A,1));
end

valBest= -inf;
val_am = -inf;
val_list=zeros(1,size(set,2));

if set==inf
    options.verbose=0;
    options.optTol=1e-8;
    options.maxIter=1000;
    options.cardinality_vec=k;
    options.initType=1; %default 2
    [uTP,val_TP] = TPower_SPCA(B, options);
    if both
        %[u ,nIter ,~, ~] = spca_am(rand(size(A,1),1), B, 2, 0, 0,k, 1000, 1e-6);
%         val_am=u'*B*u;
%         val=val_am;
        [~,Isort]=sort(diag(B),'descend');
        [vars,rhobreaks,res]=PartialPathCov(B(Isort,Isort), k);
        J=res(1:k,k);
        J=Isort(J);
        [u_am,val_am]=eigs(B(J,J),1,'lm');
        fprintf('TPower : %f am : %f\n',val_TP,val_am);        
%         keyboard;
    end
    if val_TP>val_am
        val=val_TP;
        support = find(uTP~=0);
        u=uTP(support);   
    else
        val=val_am;
        support = J;
        u=u_am;
    end
    if emin<0
        val=(val+2*emin);
    else
        val=val;
    end
    if  valBest < val
        supportBest = support;
        valBest = val;
        uBest = u;
    end
else
    for i=1:size(set,2)
        %eigenvector associated to largest real eigeinvalue
        S=find(set(:,i));
        Bs=B(S,S);
        Bs=0.5*(Bs+Bs');
        %     [v,d,flag]=eigs(-Hs,1,'la');
        [u,val,flag]=eigs(Bs,1,'lm');
        if flag
            fprintf('in get_new_atom_spca eigs has not converged.\n');
            %         d=-v'*Hs*v;
            %         keyboard;
            val=-1;
        end
        u=real(u);
        val=real(val);
        if emin<0
            val=(val+2*emin);
        else
            val=val;
        end
        val_list(i)=val;
        if val>valBest
            valBest=val;
            supportBest=S;
            uBest=u;
        end
    end
end

uBest = sparse(supportBest, ones(length(supportBest),1),uBest,size(A,1),1);
if valBest<0
    fprintf('in lmo_spca no descent direction\n');
         %keyboard;
end

