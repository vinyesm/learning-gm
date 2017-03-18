function [Dfin,Z,Z1,ActiveSet] =f2(S,lambda,mu)
    p=size(S,1);
    param.f=4;
    param.diag=0;
    param.PSD=true;
    param.max_nb_main_loop=2;%100;%2;%1000
    param.powerIter=500;
    param.stPtPowerIter=1000;
    param.niterPS=5000;%5000
    param.epsStop=1e-8;
    param.PSdualityEpsilon=1e-8;
    param.k=0;
    param.PSmu=0; %strong convexity
    param.verbose=1;
    param.debug=0;
    param.sloppy=0;
    param.max_nb_atoms=param.max_nb_main_loop*param.niterPS;
    inputData.Y=-eye(p);
    inputData.X1=S^.5;
    inputData.X2=S^.5;
    param.cardfun=inf*ones(1,p);
    param.cardfun(p)=1;
    param.lambda=lambda;
    param.mu=mu;
    [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
    if ~isempty(ActiveSet.alpha)
        Uso=bsxfun(@times,sqrt(ActiveSet.alpha)',ActiveSet.atoms);
        nl=size(ActiveSet.atoms,2);
        Dfin=zeros(p+nl);
        Dfin(1:nl,1:nl)=eye(nl);
        Dfin((nl+1):(nl+p),(nl+1):(nl+p))=-Z1;
        Dfin(1:nl,(nl+1):(nl+p))=Uso';
        Dfin((nl+1):(nl+p),1:nl)=Uso;
    else
        Dfin=Z1;
    end
end
