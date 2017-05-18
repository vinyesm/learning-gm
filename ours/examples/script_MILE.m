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

run('pp_MILE.m');

p=size(S,1);

%%



las=[.1 .05];
mus=[.001 .005 .0001];

for lam=las
    for mu=mus
        p=size(S,1);
        param.f=4;
        param.verbose=1;
        inputData.X1=S^.5;
        inputData.X2=S^.5;
        inputData.Y=-eye(p);
        

        param.cardfun=inf*ones(1,p);
        param.cardfun(100)=1;

        param.lambda=lam;
        param.mu=mu;
        param.max_nb_main_loop=100;
        [Z Z1 Z2 ActiveSet hist param flag output] = cgan_l1_omega(inputData,param);
        save(['MILE_100_lam_' num2str(lam) '_mu_' num2str(mu) '_rank_' num2str(ActiveSet.atom_count)], 'Z','Z1','Z2','ActiveSet','param');
    end
end