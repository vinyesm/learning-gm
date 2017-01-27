function [Z Znorm ActiveSet] = prox_omega(X0,lambda,options)

if nargin<3
end

p=size(X0,1);
inputData.X1=eye(p);
inputData.X2=eye(p);
inputData.Y=X0;
param.f=4;
param.PSD=true;
param.max_nb_main_loop=50;
param.powerIter=100;
param.stPtPowerIter=1000;
param.niterPS=10000;%5000
param.epsStop=1e-8;
param.PSdualityEpsilon=1e-3;
param.k=0;
param.PSmu=0; %strong convexity
param.verbose=1;
param.debug=0;
param.sloppy=0;
param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
param.cardfun=options.cardfun;
param.lambda=lambda;
param.opt='asqp';
param.diag=0;

%keyboard;
X0=(X0+X0')/2;
leig=eigs(X0,1,'la');
if leig < 0
    Z=zeros(p);
    Znorm=0;
    ActiveSet={};
else
    fprintf('cgan_spca\n');
    [Z ActiveSet hist param flag output] = cgan_spca(inputData,param);
    Znorm=sum(ActiveSet.alpha);
end

end