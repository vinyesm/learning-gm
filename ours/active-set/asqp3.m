function [c,A,nbpivot]=asqp3(Q,b,c0,param,new_atom_added,idx_atom)
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
% Resolving numerical issues with two-metric - Mark Schmidt, Dongmin Kim & Suvrit Sra
%
%
% Marina Vinyes and Guillaume Obozinski, 2016
% %%%%%%%%%%%%


max_iter=param.max_iter;
epsilon=param.epsilon;
debug_mode=param.debug_mode;
debug_mode=0;

tol=1e-12;
t=size(c0,1);
c=full(c0);
g=Q*c-b;
A=(c>tol); % active atoms
if param.ws
    if new_atom_added,
        A(idx_atom)=true; % adding the last one
        if g(idx_atom)>0,
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

iter=1;
while(iter<=max_iter)
    %% Compute new candidate solution
    J=A|(g<0);
    if norm(g(J))<epsilon,
        break;
    else
        hist.norm_g(iter)=norm(g(A));
    end
    p=zeros(t,1); 
    F=(c>tol); % free variables
    if (rcond(Q(A & F,A & F))<1e-12) %&& debug_mode,
        display('asqp : Warning, Hessian badly conditioned\n');
        keyboard;
    end
    p(A & F)=-Q(A & F,A & F)\g(A & F);
    p(A & ~F)=-g(A & ~F);
    %% Progress until active set reduces
    if (~all(p(A)+c(A)>=0)), % Drop step
%         keyboard;
        K=A & (p<0);
        [tau,i_remove]=min(c(K)./(-p(K)));
        if (tau>1),
            error('problem with tau');
        end
        I_K=find(K);
        A(I_K(i_remove))=false;
        c=c+tau*p;
        c(I_K(i_remove))=0;
        c(c<tol)=0;
        g=Q*c-b;
        A=A & (c>0);
        nb_drop_steps=nb_drop_steps+1;
%         fprintf('.');     
    else % Full step
        c=c+p;
        g=Q*c-b;
        nb_full_steps=nb_full_steps+1;
%         fprintf('+')
        if param.ws && nb_full_steps>10,
            if debug_mode
            fprintf('  nb_full_steps>10\n');
            end
            break;
        end
        if(any(g<-tol & ~A))
            [~,j]=min(g.*(~A));
            A(j)=true;
        end
    end
    if debug_mode,
        hist.obj(iter)=0.5*c'*Q*c-c'*b;
        if hist.obj(iter)>obj_old+tol,
            display('obj increases in asqp2');
%             keyboard;
        end
    end
    %% Test to increase active set
    if any(A & c==0 & g>=0)
        error('wrong feature');
    end
    if debug_mode,
        obj_old=hist.obj(iter);
        hist.J(:,iter)=A;
        hist.c(:,iter)=c;
        if(all(p(A)+c(A)>=0)),
            hist.r(iter)=0;
        else
            hist.r(iter)=I_K(i_remove);
        end
    end
    iter=iter+1;
end

nbpivot=nb_full_steps+nb_drop_steps;
%fprintf('\n');

if 0 && debug_mode,
    hist.obj=hist.obj(1:min(iter,max_iter)-1);
    hist.norm_g=hist.norm_g(1:min(iter,max_iter));
    if any(diff(hist.obj)>tol),
        display('objective increases in asqp');
    end
    if iter>max_iter || nb_full_steps>10,
        figure(15);
        plot(hist.obj);
        keyboard;
    end
end




