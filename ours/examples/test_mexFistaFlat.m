
Y=randn(10);
X=randn(10);
W0=zeros(10);
par.lambda=.1;
par.loss='cur';
par.regul='l1';
par.verbose=true;
par.max_it=10000;


[W optim]=mexFistaFlat(Y,X,W0git commit ,par);