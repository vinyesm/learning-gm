%% TEST

close all; clear all; clc;

%% fun
aff=@(a,b,x)a*x+b;

%% parameters
p=500; % vector dimention
k0=20;
u0type=1; % 1:ones, 2:randn
N=10; % nb samples for estimating expectation
sigma=.1; % noise variance
sigma2=sigma*sigma;
q=1; % BH parameter

%%
u0=[(k0:-1:1) zeros(1,p-k0)]';
% u0=[ones(1,k0) zeros(1,p-k0)]';
u0=5*u0/norm(u0);

%% Benjamini Hosberg coefficients;
kmax=p-ceil(p/2);
pkmax=p-kmax;
y=1-((1:pkmax)*q/(2*pkmax));
lBH=icdf('Normal',y,0,sigma);
l2BH=lBH.^2;

%% theoretical bound for noise
bound=cumsum(l2BH);
boundsmallk=2*sigma2*(1:p).*log(p./(1:p));

%% f candidate

% % x^pow
% f=zeros(p,1);
% ebar=l2BH(1);
% cc=4;
% slopelin=4*ebar;
% cmax=1e-16;
% pow=lambertw(cc*ebar*kmax*log(kmax/cmax))/log(kmax/cmax);
% while pow>1
%     cmax=cmax/2;
%     pow=lambertw(cc*ebar*kmax*log(kmax/cmax))/log(kmax/cmax);
% end
% f(1:kmax)=((1:kmax)./cmax).^pow;
% dkmax=f(kmax)-slopelin*kmax;
% f((kmax+1):p)=aff(slopelin,dkmax,(kmax+1):p);

% ((x-x0)/s0)^pow st f'(kmax)=ebar, f'(1)=df1
f=zeros(p,1);
ebar=l2BH(1);
%cc=4;
pow=.5;
slopelin=4*ebar;
dfkmax=4*ebar;
df1=dfkmax*2;
s0= pow^(1/pow) / ((kmax-1)^((1-pow)/pow)) * ( 1/dfkmax^(1/(1-pow)) - 1/df1^(1/(1-pow)) )^((1-pow)/pow);
x0=kmax - (pow/(dfkmax*s0^pow))^(1/(1-pow));
f(1:kmax)=(((1:kmax)-x0)./s0).^pow;
dkmax=f(kmax)-slopelin*kmax;
f((kmax+1):p)=aff(slopelin,dkmax,(kmax+1):p);

% % (x/c)^pow
% f=zeros(p,1);
% pow=.2;
% cmax=(pow/(ebar))^pow*kmax^(1-1/pow);
% f(1:kmax)=((1:kmax)./cmax).^pow;
% dkmax=f(kmax)-slopelin*kmax;
% f((kmax+1):p)=aff(slopelin,dkmax,(kmax+1):p);




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
subplot(2,2,1);
plot(Edn, 'b');hold on;
plot(Edn+Stddn, 'b--');hold on;
plot(Edn-Stddn, 'b--');hold on;
stem(k0,Edn(k0),'b--');
pbaspect([1 1 1])
subplot(2,2,2);
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
subplot(2,2,3);
plot(f, 'k-');hold on;
pbaspect([1 1 1])
subplot(2,2,4);
plot(Edn./f, 'r-');hold on;
%     stem(k0,Edn(k0)./f(k0),'b--');
stem(idx,val,'r--');
pbaspect([1 1 1])

%% Run test for all k0

%%
Ek0=zeros(p,p); %Expectation of dual norm
E2k0=zeros(p,p); %Expectation of dual norm square

%%
for k0=1:p
%     u0=[(k0:-1:1) zeros(1,p-k0)]';
    u0=[ones(1,k0) zeros(1,p-k0)]';
    u0=5*u0/norm(u0);
    for n=1:N
        x=u0+sigma*randn(p,1);
        x=sort(abs(x),'descend');
        for i=1:p
            Ek0(i,k0)=Ek0(i,k0)+sum(x(1:i).^2);
            E2k0(i,k0)=E2k0(i,k0)+sum(x(1:i).^2)^2;
        end
    end
end

xall=xall/N;
Ek0=Ek0/N;
E2k0=E2k0/N;
Vark0=(E2k0-Ek0.^2);
Stdk0=1.96*sqrt(Vark0)/sqrt(N);

[vals,idxs]=max(bsxfun(@rdivide,Ek0,f));

figure(3);
plot(1:p,idxs,'.');hold on;
plot(1:p,1:p,'r');hold on;
stem(kmax,p,'k--');
pbaspect([1 1 1]);

% figure(4);
% scatter(1:p,idxs);hold on;
% plot(1:p,1:p,'r');
% axis equal
