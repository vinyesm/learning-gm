%%
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

%%
k0=5;
p=20;
n=100;
N=100;
sigma2=.1;
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
    semilogy(1:p,Ef+Stdf,'k','LineWidth',lw); hold on;
    semilogy(1:p,Ef-Stdf,'k','LineWidth',lw); hold on;
    stem(k0,Ef(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
    
%     figure(2);clf;
%     plot(1:p,Ef,'r','LineWidth',2);
    
else
    %%
    upperEf=@(q)8*sqrt(q.*log(p./q)+2*q);
    figure(2);clf;
    plot(1:p,Ef,'LineWidth',2); hold on;
    plot(1:p,Ef+Stdf,'k','LineWidth',2); hold on;
    plot(1:p,Ef-Stdf,'k','LineWidth',2); hold on;
    %plot(1:p,upperEf(1:p),'r'); hold on;
    
end

