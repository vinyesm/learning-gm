function [cvgrid, pair_best, mincv]= cv_out(cv_cell,nb_part,pair,las,mus)

cv=zeros(nb_part,length(pair));
for j=1:nb_part
    for jj=1:length(pair)
        cv(j,jj)=cv_cell{j}{jj};
    end
end
%%
mcv=mean(cv);
cvgrid=inf*ones(length(las),length(mus));
count=1;
mincv=inf;
pair_best=0;
for i=1:length(las)
    for j=1:length(mus)
        if mus(j)>=las(i),
        cvgrid(i,j)=mcv(count);
        if mcv(count)<mincv
            mincv=mcv(count);
            pair_best=count;
        end
        count=count+1;
        end
    end
end