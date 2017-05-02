function [ I, beta ] = mat2l1index(M,l1_atoms)
% for symmetric matrices

p=size(M,1);
I=[];
beta=[];

for row=1:p
    for col=1:p
        val=M(row,col);
        if val~=0
            sa=sign(val);
            i1=(row-1)*p+col;
            i2=(col-1)*p+row;
            idx_l1 = find((l1_atoms(i1, :) == sa) & (l1_atoms(i2, :) == sa));
            I=[I idx_l1];
            beta=[beta ; abs(val)];
        end
    end
end

end

