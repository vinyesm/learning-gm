clear all;
addpath('../main');
addpath('../active-set');
addpath('../atom-selection');
addpath('../utils');
% addpath('../other');
addpath('../prox');
addpath('../TPower_1.0');
addpath('../TPower_1.0/algorithms/TPower/');
addpath('../TPower_1.0/misc/');
load('Lsl.mat');


las=[1 .1 .01];
ks=[200 150 100 50 inf];

for lam=las
    for k=ks
        p=size(Lsl,1);
        param.f=5;
        param.verbose=1;
        inputData.X=inv(2*eigs(Lsl,1, 'la')*eye(p)-Lsl);
        inputData.Y=-eye(p);
        
        if k==inf
            beta=.5;
            param.cardfun=((1:p).^beta)/p^beta;
            param.cardfun(1)=inf;
            param.cardfun(200:end)=inf;
        else
            param.cardfun=inf*ones(1,p);
            param.cardfun(k)=1;
        end
        param.lambda=lam;
        param.mu=.01;
        param.max_nb_main_loop=100;
        [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
        save(['k_' num2str(k) '_lam_' num2str(lam) '_rank_' num2str(ActiveSet.atom_count)], 'Z','Z1','Z2','ActiveSet','param');
    end
end