close all; clear all; clc;

%%
p=100;
k0=20;
dn=k0;
a1=1; % 0<a1<=1
a2=0.5; % 0<a2<a1

%%
aff=@(a,b,x)a*x+b;

%% we want a1*k0+b1=a2*k0+b2 and a1*k0+b1=dn
b1=dn-a1*k0;
b2=dn-a2*k0;
x1=1:k0;
x2=(k0+1):p;
h=[aff(a1,b1,x1) aff(a2,b2,x2)];

%% candidat ex+d
e=.1;
d=e/2*(b1/a1+b2/a2);
f=aff(e,d,1:p);

%%
figure(1);
subplot(1,2,1);
plot(h,'b'); hold on;
plot(f,'g'); hold on;
subplot(1,2,2);
plot(h./f,'r'); hold on;