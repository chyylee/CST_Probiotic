function [CST_post] = get_CST_post(yout)

%%
    SSref = [0.911801263	0.059240659	0.028958078;
        0.146456352	0.758938299	0.09460535;
        0.153186234	0.09522381	0.751589956];
    
    rel_out = yout ./ sum(yout,2);
    
    run_mat = [rel_out(:,1:2) rel_out(:,3)+rel_out(:,4)];
    
    for i = 1:size(run_mat,1)
        if sum(run_mat(i,:) < -1E-3) > 0 % if values are less than zero
            CST_post(i) = NaN;
        else
    
            dif = sum((SSref - repmat(run_mat(i,:),size(SSref,1),1)).^2,2);
            [~,id] = min(dif);
            
            if id == 3
                if yout(i,3) <= yout(i,4)
                    id = 4;
                else
                    id = 3;
                end
            end
            CST_post(i) = id;
        end
    end
end