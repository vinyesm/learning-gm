clear all; clc;

load('dualnorm-noise-ub.mat');

upperEf4=@(q)4*sqrt(q.*log(p./q)+2*q); %tighter for symmetric matrices
lw=1;

figure(1);clf;
plot(1:p,Ef,'r','LineWidth',lw); hold on;
jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
xlabel('k');
ylabel('E[\Omega_k^*(XX'')]')
title('Expectation of dual norm of noise XX'' for X standard iid ')
legend(['E[\Omega_k^*(XX'')]'], 'Location','NorthWest');
set(gca,'fontsize',15);
print('-depsc','-tiff','dualnorm-noise')

figure(2);clf;
plot(1:p,Ef,'r','LineWidth',lw); hold on;
jbfill(1:p,Ef+Stdf,Ef-Stdf,ones(p,1),'r','r',1,.1);hold on;
plot(1:p,upperEf4(1:p),'Color','k','LineWidth',lw); hold on;
xlabel('k');
ylabel('E[\Omega_k^*(XX'')]')
title('Expectation of dual norm of noise XX'' for X standard iid ')
legend({'E[\Omega_k^*(XX'')]', 'variance','upper bound'}, 'Location','NorthWest');
set(gca,'fontsize',15);
print('-depsc','-tiff','dualnorm-noise-ub')