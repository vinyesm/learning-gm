function [c,A,nbpivot,ng]=asqpnb(Q,b,c0,param,new_atom_added,idx_atom,atom)
%
% SOLVING QUADRATIC PROBLEM WITH ACTIVE SET, WITH WARM START
% min_c 0.5*c'Qc-c'*b s.t. c \geq 0

% INPUTS:
% Q,b,c0           = quadratic problem parameters and intitial point c0
% new_atom_added   = weather or not a new atom has been added to the active set
% param.max_iter   = max nb of iterations
% param.epsilon    = tolerance parameter
% param.debug_mode = true/false
%
% OPTIONAL INPUTS
% param.max_nb_iter    = default 500 = maximum number of iterations allowed default = 500;
% param.max_nb_atoms   = default 500 = maximum number of atoms (used to allocate memory)
% param.epsStop        = default 1e-5 = tolerance parameter
% param.flag_no_design = when no design matrix more optimized storage.
% param.debug          = false = to debug;
% param.debug_asqp     = false = to debug inside asqp
%
% OUTPUTS :
% c       = final result
% A       = vector of 0/1 indicating active atoms
% nbpivot = final nb of pivots, i.e. full_steps+drop_steps
%
%%%%%%%%%%%%
%
% Active-set algorithm from
% Numerical Optimization - J.Nocedal & S.J.Wright
%
% Marina Vinyes and Guillaume Obozinski, 2016
% %%%%%%%%%%%%

% keyboard

MAX_NB_FULL_STEPS=200;
max_iter=param.max_iter;
epsilon=param.epsilon;
debug_mode=param.debug_mode;
debug_mode=0;

tol=1e-12;
epsilon=tol*10;
t=size(c0,1);
c=full(c0);
g=Q*c-b;
A=(c>tol); % active atoms
if param.ws
    if new_atom_added,
        %         A(idx_atom)=true; % adding the last one
        j=idx_atom;
        if g(j)>0,
            error('The new atom direction is not a descent direction')
        end
    end
else
    if(any(g<-tol))
        [~,j]=min(g);
        A(j)=true;
    else
        error('No direction added')
    end
end

hist.norm_g=zeros(1,max_iter);
nb_drop_steps=0;
nb_full_steps=0;

if debug_mode,
    hist.obj=zeros(1,max_iter);
    hist.J=zeros(t,max_iter);
    hist.c=zeros(t,max_iter);
    hist.r=zeros(1,max_iter);
    obj_old=0.5*c'*Q*c-c'*b;
end

E=eye(t);
iter=1;
while(iter<=max_iter)
    if g(j)>0
        break;
    else
        if A(j)==true
            error('error j should not belong to the working set');
        end
        C=~A;
        C(j)=true;
        jj=find(C==j);
        m=sum(C);
        KKT = zeros(t+m,t+m);
        KKT(1:t,1:t)=Q;
        KKT(1:t,t+1:end)=E(:,C);
        KKT(t+1:end,1:t)=E(C,:);
        r=KKT\[zeros(t,1);E(C,j)];
        p=r(1:t);
        q=-r(t+1:end);
        K=A & (p<0);
        [tauF,i_remove]=min(c(K)./(-p(K)));
        I_K=find(K);
        if (tauF>1),
            error('problem with tau');
        end
        tauS=-g(j)/q(jj);
        if tauS<0
            error('unbounded');
        end
        tau=min(tauF,tauS);
        c=c+tau*p;
        g=g*tau*Q*p;
        if tauF<tauS
            % remove constraint, increase active set
            r=KKT\[E(:,I_K(i_remove));zeros(m,1)];
            u=r(1:t);
            v=-r(t+1:end);
            if norm(u)<1e-16
                sig=-g(j)/q(jj);
            else
                sig=0;
            end
            g
        else
            %
        end
        iter=iter+1;
    end
end

nbpivot=nb_full_steps+nb_drop_steps;
ng=norm(g(J));
%fprintf('\n');

if debug_mode %&& 0,
    hist.obj=hist.obj(1:min(iter,max_iter)-1);
    hist.norm_g=hist.norm_g(1:min(iter,max_iter));
    if any(diff(hist.obj)>tol),
        display('objective increases in asqp');
    end
    if iter>max_iter || nb_full_steps>MAX_NB_FULL_STEPS,
        figure(15);
        plot(hist.obj);
        keyboard;
    end
end




