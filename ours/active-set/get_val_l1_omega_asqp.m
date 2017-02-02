function [loss,pen,obj,dualityGap,time]=get_val_l1_omega_asqp(ActiveSet,inputData,param,cardVal)

currloss = .5*norm(inputData.Y - inputData.X1*Z*inputData.X2, 'fro')^2;
loss= currloss;
pen= param.lambda*dot(ActiveSet.alpha,cardVal)+param.mu*sum(ActiveSet.beta);
obj = currloss + pen ;
dualityGap = get_dg_l1_omega_asqp(Z,param,ActiveSet,inputData,obj,loss,pen);
dualityGap= real(dualityGap);
time=toc;
inputData.Y= YStart;

end

