%% TEST

close all; clear all; clc;

%% fun
aff=@(a,b,x)a*x+b;

%% parameters
p=100; % vector dimention
k0=20;
u0type=1; % 1:ones, 2:randn
N=10; % nb samples for estimating expectation
sigma=.1; % noise variance
sigma2=sigma*sigma;
q=1; % BH parameter

%%
u0=[(k0:-1:1) zeros(1,p-k0)]';
u0=5*u0/norm(u0);

%% Benjamini Hosberg coefficients;
kmax=p-ceil(p/3);
pkmax=p-kmax;
y=1-((1:pkmax)*q/(2*pkmax));
lBH=icdf('Normal',y,0,sigma);
l2BH=lBH.^2;

%% theoretical bound for noise
bound=cumsum(l2BH);
boundsmallk=2*sigma2*(1:p).*log(p./(1:p));

%% f candidate

% x^pow
ebar=l2BH(1);
cc=4;
slopelin=4*ebar;
cmax=1e-16;
pow=lambertw(cc*ebar*kmax*log(kmax/cmax))/log(kmax/cmax);
while pow>1
    cmax=cmax/2;
    pow=lambertw(cc*ebar*kmax*log(kmax/cmax))/log(kmax/cmax);
end


% % (x/c)^pow
% pow=.2;
% cmax=(pow/(ebar))^pow*kmax^(1-1/pow);

f=zeros(p,1);
f(1:kmax)=((1:kmax)./cmax).^pow;
% f(1:kmax)=log(1+(1:kmax));
dkmax=f(kmax)-slopelin*kmax;
f((kmax+1):p)=aff(slopelin,dkmax,(kmax+1):p);


%%
Edn=zeros(p,1); %Expectation of dual norm
Edn2=zeros(p,1); %Expectation of dual norm square
xall=zeros(p,1);

%% 
for n=1:N
    x=u0+sigma*randn(p,1);
    xall=xall+x;
    x=sort(abs(x),'descend');
    for i=1:p
        Edn(i)=Edn(i)+sum(x(1:i).^2);
        Edn2(i)=Edn2(i)+sum(x(1:i).^2)^2;
    end
end

xall=xall/N;
Edn=Edn/N;
Edn2=Edn2/N;
Vardn=(Edn2-Edn.^2);
Stddn=1.96*sqrt(Vardn)/sqrt(N);

[val,idx]=max(Edn./f);

figure(1);
bar(abs(xall));

figure(2);
subplot(1,3,1);
    plot(Edn, 'b');hold on;
    plot(Edn+Stddn, 'b--');hold on;
    plot(Edn-Stddn, 'b--');hold on;
    %plot(bound,'r'); hold on;
    %plot(boundsmallk,'g'); hold on;
    plot(f, 'k-');hold on;
    plot(Edn./f, 'r-');hold on;
    stem(k0,Edn(k0),'b--');
    stem(idx,val,'r--');
    pbaspect([1 1 1])
subplot(1,3,2);
    plot(f, 'k-');hold on;
    pbaspect([1 1 1])
subplot(1,3,3);
    plot(Edn./f, 'r-');hold on;
%     stem(k0,Edn(k0)./f(k0),'b--');
    stem(idx,val,'r--');
    pbaspect([1 1 1])
