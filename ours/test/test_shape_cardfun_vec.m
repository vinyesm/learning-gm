%%
close all; clear all; clc;

%% parameters
p=100; % vector dimention
N=10; % nb samples for estimating expectation
sigma2=1; % noise variance
q=1; % BH parameter


%% Benjamini Hosberg coefficients;
y=1-(1:p)*q/(2*p);
lBH=icdf('Normal',y);
l2BH=lBH.^2;

%% theoretical bound for noise
bound=cumsum(l2BH);

%%
Edn=zeros(p,1); %Expectation of dual norm
Edn2=zeros(p,1); %Expectation of dual norm square

%% 
for n=1:N
    x=randn(p,1);
    x=sort(abs(x),'descend');
    for i=1:p
        Edn(i)=Edn(i)+sum(x(1:i).^2);
        Edn2(i)=Edn2(i)+sum(x(1:i).^2)^2;
    end
end

Edn=Edn/N;
Edn2=Edn2/N;
Vardn=(Edn2-Edn.^2);
Stddn=1.96*sqrt(Vardn)/sqrt(N);

figure(1);
plot(Edn, 'b');hold on;
plot(Edn+Stddn, 'b--');hold on;
plot(Edn-Stddn, 'b--');hold on;
plot(bound,'r'); hold on;