clear all; clc;

k0=5;
% savename=['dualnorm-u0rand-k0-' num2str(k0)];
savename=['dualnorm-u0ones-k0-' num2str(k0)];
load([savename '.mat']);

%% upper bounds
upperEf=@(q)8*sqrt(q.*log(p./q)+2*q);
upperEf2=@(q)32*(q.*log(p./q)+2*q);
upperEf4=@(q)4*sqrt(q.*log(p./q)+2*q); %tighter for symmetric matrices
upper3=@(q,t)sqrt(log(nchoosek(p,q))/t+16*q+8/(1-8*t));
upperEf3=zeros(p,3);
for q=1:p
    upperEf3(q,1)=upper3(q,.12);
end;

%% candidate fun
% param for candidate fun
aaa=.1;
max2u0=max(u0.^2);
max2u0=max2u0(1);
norm2u0=norm(u0)^2;
a0=(norm2u0-max2u0)/(k0^aaa-1);
b0=max2u0-a0;
ap=(norm2u0-max2u0)/((p/2)^aaa-1);
bp=max2u0-a0;
crand=a0;
drand=max(0,a0*a0/(b0+norm2u0));
candidate2=@(x)crand*x.^aaa+drand;
% param for candidate fun
candidate=@(x)norm(u0)^2*(2*x/p+1);

%% plots
lw=1;
figure(1);clf;
plot(1:p,Ef,'r','LineWidth',lw); hold on;
jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
stem(k0,Ef(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
xlabel('k');
ylabel('E[\Omega_k^*(M)]')
title('Expectation of dual norm of M')
legend({'Expectation of dual norm'}, 'Location','NorthWest');
set(gca,'fontsize',15);
print('-depsc','-tiff',savename)

figure(2);clf;
plot(1:p,Ef,'r','LineWidth',lw); hold on;
jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
stem(k0,Ef(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
plot(1:p,upperEf3(:,1),'Color',[0,1,1],'LineWidth',lw); hold on;
plot(1:p,upperEf4(1:p),'Color','k','LineWidth',lw); hold on;
plot(1:p,Ef-upperEf3(:,1)','LineStyle','-.','Color',[0,1,1],'LineWidth',lw); hold on;
plot(1:p,Ef-upperEf4(1:p),'LineStyle','-.','Color','k','LineWidth',lw); hold on;
xlabel('k');
%ylabel('E[\Omega_k^*(\hat{M})]')
title('Expectation of dual norm of M')
legend({'Expectation of dual norm' 'avriance', 'ub1','ub2', 'E-ub1','E-ub2'}, 'Location','SE');
set(gca,'fontsize',15);
print('-depsc','-tiff',[savename '-minusnoise'])


[v1, i1]=max(Ef./candidate(1:p));
[v2, i2]=max(Ef./candidate2(1:p));

figure(4);clf;
plot(1:p,Ef./candidate(1:p),'Color',[0,1,0.2],'LineWidth',lw); hold on;
jbfill(1:p,(Ef+Stdf)./candidate(1:p),(Ef-Stdf)./candidate(1:p),ones(p,1),'g','g',1,.1);hold on;
stem(k0,Ef(k0)/candidate(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
stem(i1,v1,'LineStyle','-.','Color',[0,1,0],'LineWidth',lw);hold on;
pbaspect([1 1 1]);
xlabel('k');
ylabel('Expectation of dual norm/candidate')
title('linear candidate')
legend({'E/candidate' 'variance'}, 'Location','SE');
set(gca,'fontsize',15);
print('-depsc','-tiff',[savename '-candidatelin'])

figure(4);clf;
plot(1:p,Ef./candidate2(1:p),'Color',[0,1,0.2],'LineWidth',lw); hold on;
jbfill(1:p,(Ef+Stdf)./candidate2(1:p),(Ef-Stdf)./candidate2(1:p),ones(p,1),'g','g',1,.1);hold on;
stem(k0,Ef(k0)/candidate2(k0),'LineStyle','-.','Color',[1,0,0],'LineWidth',lw);hold on;
stem(i2,v2,'LineStyle','-.','Color',[0,1,0],'LineWidth',lw);hold on;
pbaspect([1 1 1]);
xlabel('k');
ylabel('Expectation of dual norm/candidate')
title(['pow ' num2str(aaa) ' candidate'])
legend({'E/candidate' 'variance'}, 'Location','SE');
set(gca,'fontsize',15);
print('-depsc','-tiff',[savename '-candidatepow'])
