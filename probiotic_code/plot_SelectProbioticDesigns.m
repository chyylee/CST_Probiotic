function [allCounts,pChiDfAllResponses, pChiDfBvResponses] = plot_SelectProbioticDesigns(pidx_list,...
    pert_list,pert_names,gr1,gr2,all_select_outcomes,all_warnall,time_id)

    [param_names,~,sp_cols,~] = get_naming_terms();
    
    % PULL DATA FOR PLOTTING 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    allCounts = NaN(length(pidx_list),4);
    allFrequencies = NaN(length(pidx_list),4);
    param_vals = zeros(7,length(pidx_list));
    for pert_id = 1:length(pidx_list)
        abundance = squeeze(all_select_outcomes(pert_id,all_warnall(pert_id,:)==1,time_id,2:5));
        CST_post = get_CST_post(abundance);
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts,2);
        allCounts(pert_id,:) = Counts;
        allFrequencies(pert_id,:) = Frequency;
        
        param_vals(pidx_list{pert_id},pert_id) = pert_list{pert_id};
    end
    
   
    % STATISTICS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pChiDfAllResponses = NaN(length(pidx_list),length(pidx_list),3);
    pChiDfBvResponses = NaN(length(pidx_list),length(pidx_list),3);
    for pert_id1 = 1:length(pidx_list)
        for pert_id2 = 1:length(pidx_list)
            observed = [allCounts(pert_id1,:)', allCounts(pert_id2,:)'];
            [p, chi2stat,df] = chigof(observed);
            pChiDfAllResponses(pert_id1,pert_id2,:) = [p, chi2stat,df]; % Across all 4 reseponse types
        
            observed_BV = [observed(1,:);sum(observed(2:end,:))];
            [p_BV, chi2stat_BV,df_BV] = chigof(observed_BV);
            pChiDfBvResponses(pert_id1,pert_id2,:) = [p_BV, chi2stat_BV,df_BV]; % BV+/-
        end
    end
    
    
    % PLOT THE RESULTS 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ~~~~~~~~~~~~ Plots Parameter Alterations  ~~~~~~~~~~~~ 
    subplot(4,2,7)
    param_vals(param_vals > 1) = 1;
    param_vals(isnan(param_vals)) = 0;
    ax = imagesc(param_vals([1,2,5],gr1));
    colormap(gca,[0 0.75 0.75; ...
        1 1 1; ...
        0.75 0 0.75])
    set(gca,'YTICK',1:3,'YTICKLABELS',param_names([4,8,17]),...
        'XTICK',1:length(gr1),'XTICKLABELS',pert_names(gr1))
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    xtickangle(90)
    
    subplot(4,2,8)
    param_vals(param_vals > 1) = 1;
    ax = imagesc(param_vals([5,6,7],gr2));
    colormap(gca,[0 0.75 0.75; ...
        1 1 1; ...
        0.75 0 0.75])
    set(gca,'YTICK',1:3,'YTICKLABELS',param_names([17,18,19]),...
        'XTICK',1:length(gr2),'XTICKLABELS',pert_names(gr2))
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    xtickangle(90)
    
    % ~~~~~~~~~~~~ Plots the BV Design Comparisons  ~~~~~~~~~~~~ 
    subplot(4,2,[1,3,5])
    X = allFrequencies(gr1,:)*100;
    bar(X,'stacked')
    colororder(sp_cols)
    xlim([0.5 length(gr1)+0.5])
    xticklabels([])
    ylabel('Percent Samples')
    hold on
    Y = ones(length(gr1),1)*100;
    P = squeeze(pChiDfBvResponses(gr1,gr1,1));
    superbar(Y, 'P',round(P,4),'BarFaceColor', 'none', 'BarEdgeColor', 'none')
    ylim([0 220])
    
    for i = 1:size(X,1) % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f %%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center',...
            'fontsize',10)
    end
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    title('Reducing BV+')
    
    % ~~~~~~~~~~~~ Plots the LB Design Comparisons  ~~~~~~~~~~~~ 
    subplot(4,2,[2,4,6])
    X = allFrequencies(gr2,:)*100;
    bar(X,'stacked')
    colororder(sp_cols)
    xlim([0.5 length(gr2)+0.5])
    xticklabels([])
    ylabel('Percent Samples')
    hold on
    Y = ones(length(gr2),1)*100;
    P = squeeze(pChiDfBvResponses(gr2,gr2,1));
    superbar(Y, 'P',round(P,4),'BarFaceColor', 'none', 'BarEdgeColor', 'none')
    ylim([0 220])
    
    for i = 1:size(X,1) % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f %%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center',...
            'fontsize',10)
    end
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    title('Boosting oLB')
end