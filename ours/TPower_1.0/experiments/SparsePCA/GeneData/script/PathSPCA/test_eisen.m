% Test code on various examples
tic;clear all;



load ..\..\data\EisenData.mat
S=eisen2+eye(size(eisen2))*1e-14;
[d,ix]=sort(diag(S),'descend');S=S(ix,ix);
F=chol(S);
k=size(S,1);
kp=k;

[vars,rhobreaks,res]=FullPathCov(S);
t_pathspca = toc
[bnds,rhov,dualvals]=UpperBounds(F,S,res(:,1:kp));


% Plot with optimal points. Uncomment to merge multiple plots
epsr=1e-4;
error=max(bnds(1:kp)-vars(1:kp),zeros(size(vars(1:kp))));
optpi=find((error./vars(1:kp))<=epsr);optpv=vars(find((error./vars(1:kp))<=epsr));
%errorl=max(bndsl-varsl(1:kp),zeros(size(vars(1:kp))));
%optpil=find((errorl./varsl(1:kp))<=epsr);optpvl=varsl(find((errorl./varsl(1:kp))<=epsr));
plot(1:kp,vars(1:kp),'-b',1:kp,bnds,':k','LineWidth',1.5);hold on;
%plot(1:kp,varsl(1:kp),'-b',1:kp,bndsl,':k','LineWidth',1.5);
%plot(optpil,optpvl,'b.','MarkerSize',20);
plot(optpi,optpv,'b.','MarkerSize',20);
xlabel('card');ylabel('var');
%semilogy(error./vars,'bo');xlabel('Card');ylabel('Relative Error');axis([1 kp 1e-6 1e2]);

variance_rate = vars / trace(eisen2);
bnds_rate = bnds / trace(eisen2);

% save results
PathSPCA_Result_eisen = [];
PathSPCA_Result_eisen.CPUTime = t_pathspca;
PathSPCA_Result_eisen.variance_rate = variance_rate;
PathSPCA_Result_eisen.variance_rate_bound = bnds_rate;
PathSPCA_Result_eisen.cardinality_path = 1:kp;
PathSPCA_Result_eisen.optpi = optpi;

save ..\..\result\PathSPCA_Result_eisen.mat 'PathSPCA_Result_eisen';