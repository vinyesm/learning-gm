function []=check_complementary_slackness(init,set,val_list,param)
for i=1:size(set,2)
    atomidxset_for_i=find(all(bsxfun(@eq,set(:,i),(init.atoms_u(:,1:init.atomCount)~=0))));
    if ~isempty(atomidxset_for_i)
        Z=init.atoms_u(:,atomidxset_for_i)*diag(init.coeff(atomidxset_for_i))*init.atoms_u(:,atomidxset_for_i)';
        if abs(norm(full(Z))*(param.lambda-val_list))./param.lambda>1e-1
            fprintf('complementary slackness does not hold\n');
            %keyboard
        end
    end
end