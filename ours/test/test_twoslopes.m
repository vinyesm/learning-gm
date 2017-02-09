close all; clear all; clc;

%%
p=100;
k0=20;
dn=k0;
a1=.5; % 0<a1<=1
a2=0.1; % 0<a2<a1

%%
aff=@(a,b,x)a*x+b;

%% we want a1*k0+b1=a2*k0+b2 and a1*k0+b1=dn
b1=dn-a1*k0;
b2=dn-a2*k0;
x1=1:k0;
x2=(k0+1):p;
h=[aff(a1,b1,x1) aff(a2,b2,x2)];

%% candidat ex+d
%e=1;
% d=e/2*(b1/a1+b2/a2);
e=.5*(a1+a2);
d=1;
% b1c=1;
% b2c=1;
% d=e/2*(b1c/a1+b2c/a2);
f=aff(e,d,1:p);

%%

fprintf('left k0 derivative of h/f %f\n',a1*d-e*b1); %(a1*d-e*b1) > 0
fprintf('right k0 derivative of h/f %f\n',a2*d-e*b2); %(a2*d-e*b2) < 0
fprintf('lb for coeff d : %f\n' ,e*b1/a1)
fprintf('ub for coeff d : %f\n' ,e*b2/a2)

%%
[mhf,ihf]=max(h./f(1:p));

figure(1);
subplot(1,2,1);
plot(h,'b'); hold on;
stem(k0,h(k0),'b--');
plot(f,'g'); hold on;
subplot(1,2,2);
plot(h./f,'r'); hold on;
stem(ihf,mhf,'r--');