%%
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

%%
k0=0; % 0 for only noise
p=10;
n=50;
N=10;
sigma2=1;
u0=[ones(k0,1);zeros(p-k0,1)];
u0=u0/norm(u0);

if k0>0
    M=u0*u0'+sigma2*eye(p);% covariance with noise
else
    M=sigma2*eye(p); % only noise
end



%%
%X=mvnrnd(zeros(p,1), M, n)';
%S=cov(X');

Ef=zeros(1,p);
Ef2=zeros(1,p);

for i=1:N
    fprintf('N=%d/%d\n',i,N);
    X=mvnrnd(zeros(p,1), M, n)';
    S0=sigma2*eye(p);
    S=X*X'/n;
    for k=1:p
        options.verbose=0;
        options.optTol=1e-8;
        options.maxIter=1000;
        options.cardinality_vec=k;
        [u,omegastar] = TPower_SPCA(S, options);
        Ef(k)=Ef(k)+omegastar;
        Ef2(k)=Ef2(k)+omegastar^2;
    end
end

Ef=Ef/N;
Ef2=Ef2/N;
Varf=(Ef2-Ef.^2);
Stdf=1.96*sqrt(Varf)/sqrt(N);

if k0>0
    %%
    lw=1;
    figure(1);clf;
%     errorbar(1:p,Ef,Stdf,Stdf)
    semilogy(1:p,Ef,'r','LineWidth',lw); hold on;
    jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
    stem(k0,Ef(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
    
%     figure(2);clf;
%     plot(1:p,Ef,'r','LineWidth',2);
    
else
    %%
    upperEf=@(q)8*sqrt(q.*log(p./q)+2*q);
    upperEf2=@(q)32*(q.*log(p./q)+2*q);
    upper3=@(q,t)sqrt(log(nchoosek(p,q))/t+16*k+8/(1-8*t));
    
    upperEf3=zeros(p,3);
    for q=1:p
        upperEf3(q,1)=upper3(q,.12);
        upperEf3(q,2)=upper3(q,.1);
        upperEf3(q,3)=upper3(q,10);
    end;
    
    figure(1);clf;
    semilogy(1:p,Ef,'r','LineWidth',2); hold on;
    jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
    
    figure(2);clf;
    semilogy(1:p,Ef,'r','LineWidth',2); hold on;
    jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
    %semilogy(1:p,upperEf(1:p),'k'); hold on;
    %semilogy(1:p,upperEf2(1:p),'b'); hold on;
    semilogy(1:p,upperEf3(:,1),'Color',[0,0,1],'LineWidth',2); hold on;
    semilogy(1:p,upperEf3(:,2),'Color',[0,.6,1],'LineWidth',2); hold on;
    semilogy(1:p,upperEf3(:,3),'Color',[0,1,1],'LineWidth',2); hold on;
    
end

