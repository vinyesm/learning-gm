close all; clear all; clc;

%%
p=100;
k0=80;
dn=k0;
a1=.2; % 0<a1<=1
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
%
% comment : maybe only condition is a2/b2 < e/d < a1/b1
%
% th choice (knowing k0)

%% affine candidate knowing k0
%e=1;
% d=e/2*(b1/a1+b2/a2); % (*)
% d=e/2*((dn/a1+dn/a2)-2*k0); % equivalent to (*)
% d=e/2*(b1/a1+b1/a2 +(a1/a2-1)*k0); % equivalent to (*)

%% affine candidate, independent of k0
e=1; % any >0
b1=100;
d=e/2*(b1/a1+b1/a2 +(a1/a2-1)*1); % k0=1, candidate

%%

fprintf('left k0 derivative of h/f %f\n',a1*d-e*b1); %(a1*d-e*b1) > 0
fprintf('right k0 derivative of h/f %f\n',a2*d-e*b2); %(a2*d-e*b2) < 0
fprintf('lb for coeff d : %f\n' ,e*b1/a1)
fprintf('ub for coeff d : %f\n' ,e*b2/a2)

%%
f=aff(e,d,1:p);
[mhf,ihf]=max(h./f(1:p));

figure(1);
subplot(1,2,1);
plot(h,'b'); hold on;
stem(k0,h(k0),'b--');
plot(f,'g'); hold on;
subplot(1,2,2);
plot(h./f,'r'); hold on;
stem(ihf,mhf,'r--');