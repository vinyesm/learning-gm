%%
clear all; clc;
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');

%%
RAND_U0=true;

%%
k0=80; % 0 for only noise
p=100;
n=1000;
N=40;
sigma2=1;

if RAND_U0
    u0=[rand(k0,1);zeros(p-k0,1)];
else
    u0=[ones(k0,1);zeros(p-k0,1)];
end
u0=5*u0/norm(u0);

if k0>0
    M=u0*u0'+sigma2*eye(p);% covariance with noise
else
    M=sigma2*eye(p); % only noise
end

%%candidate fun
% param for candidate fun
aaa=.1;
max2u0=max(u0.^2);
max2u0=max2u0(1);
norm2u0=norm(u0)^2;
a0=(norm2u0-max2u0)/(k0^aaa-1);
b0=max2u0-a0;
ap=(norm2u0-max2u0)/((p/2)^aaa-1);
bp=max2u0-a0;
% crand=a0;
% drand=max(0,a0*a0/(b0+norm2u0));
% crand=1/norm2u0;
% drand=crand*b0/a0+norm2u0*crand;
% crand=ap;
% drand=max(0,ap*ap/(bp+norm2u0));
crand=1/norm2u0;
drand=crand*bp/ap+norm2u0*crand;
candidate2=@(x)crand*x.^aaa+drand;
% param for candidate fun
candidate=@(x)norm(u0)^2*(2*x/p+1);

%% upper bounds
upperEf=@(q)8*sqrt(q.*log(p./q)+2*q);
upperEf2=@(q)32*(q.*log(p./q)+2*q);
upperEf4=@(q)4*sqrt(q.*log(p./q)+2*q); %tighter for symmetric matrices
upper3=@(q,t)sqrt(log(nchoosek(p,q))/t+16*q+8/(1-8*t));
upperEf3=zeros(p,3);
for q=1:p
    upperEf3(q,1)=upper3(q,.12);
    upperEf3(q,2)=upper3(q,.1);
    upperEf3(q,3)=upper3(q,.01);
end;



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
lw=1;
if k0>0
    %%
    figure(1);clf;
    plot(1:p,Ef,'r','LineWidth',lw); hold on;
    jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
    stem(k0,Ef(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
    
    figure(2);clf;
    plot(1:p,Ef,'r','LineWidth',lw); hold on;
    jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
    stem(k0,Ef(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
    plot(1:p,upperEf3(:,1),'Color',[0,1,1],'LineWidth',lw); hold on;
    plot(1:p,upperEf4(1:p),'Color','k','LineWidth',lw); hold on;
    plot(1:p,Ef-upperEf3(:,1)','LineStyle','-.','Color',[0,1,1],'LineWidth',lw); hold on;
    plot(1:p,Ef-upperEf4(1:p),'LineStyle','-.','Color','k','LineWidth',lw); hold on;
    
    figure(3);clf;
    [v1, i1]=max(Ef./candidate(1:p));
    [v2, i2]=max(Ef./candidate2(1:p));
    subplot(1,2,1)
    plot(1:p,Ef./candidate(1:p),'Color',[0,1,0.2],'LineWidth',lw); hold on;
    jbfill(1:p,(Ef+Stdf)./candidate(1:p),(Ef-Stdf)./candidate(1:p),ones(p,1),'g','g',1,.1);hold on;
    stem(k0,Ef(k0)/candidate(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
    stem(i1,v1,'LineStyle','-.','Color',[0,1,0],'LineWidth',lw);hold on;
    pbaspect([1 1 1]);
    subplot(1,2,2)
    plot(1:p,Ef./candidate2(1:p),'Color',[0,1,0.2],'LineWidth',lw); hold on;
    jbfill(1:p,(Ef+Stdf)./candidate2(1:p),(Ef-Stdf)./candidate2(1:p),ones(p,1),'g','g',1,.1);hold on;
    stem(k0,Ef(k0)/candidate2(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
    stem(i2,v2,'LineStyle','-.','Color',[0,1,0],'LineWidth',lw);hold on;
    pbaspect([1 1 1]);
    
    if RAND_U0
        save(['dualnorm-u0rand-k0-' num2str(k0)],'k0','p','n','N','sigma2','u0','M','Ef','Ef2','Varf','Stdf');
    else
        save(['dualnorm-u0ones-k0-' num2str(k0)],'k0','p','n','N','sigma2','u0','M','Ef','Ef2','Varf','Stdf');
    end
else
    
    figure(1);clf;
    plot(1:p,Ef,'r','LineWidth',lw); hold on;
    jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
    
    figure(2);clf;
    plot(1:p,Ef,'r','LineWidth',lw); hold on;
    jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
    %     plot(1:p,upperEf3(:,1),'Color',[0,0,1],'LineWidth',lw); hold on;
    %     plot(1:p,upperEf3(:,2),'Color',[0,.6,1],'LineWidth',lw); hold on;
    %     plot(1:p,upperEf3(:,3),'Color',[0,1,1],'LineWidth',lw); hold on;
    plot(1:p,upperEf4(1:p),'Color','k','LineWidth',lw); hold on;
    
%     save('dualnorm-noise-ub','k0','p','n','N','sigma2','u0','M','Ef','Ef2','Varf','Stdf')
    
end

