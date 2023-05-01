function [allCounts,pChiDfAllResponses, pChiDfBvResponses] = ...
    plot_4D_SensitivityAnalysis(all_select_outcomes,all_warnall,...
    newValueMat,pidx,time_id,param_names,sp_cols)
    % PULL DATA FOR PLOTTING 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    allCounts = NaN(size(newValueMat,1),4);
    allFrequencies = NaN(size(newValueMat,1),4);
    for pert_id = 1:size(newValueMat,1)
        Xwarn = all_warnall(pert_id,:) == 1;
        X = squeeze(all_select_outcomes(pert_id,Xwarn,time_id,2:end));
        CST_post = get_CST_post(X);
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts,2);
        allCounts(pert_id,:) = Counts;
        allFrequencies(pert_id,:) = Frequency;
    end
    
    [v,idx_sorted] = sort(allFrequencies(:,1),'ascend'); % Best to worst in nAB response frequency
    
    % STATISTICS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    base_id = find(sum(abs(newValueMat),2) == 0); % Baseline vs all others
    base_dat = allCounts(base_id,:)';
    
    pChiDfAllResponses = NaN(size(newValueMat,1),3);
    pChiDfBvResponses = NaN(size(newValueMat,1),3);
    
    for pert_id = 1:size(newValueMat,1)
        observed = [base_dat, allCounts(pert_id,:)'];
        [p, chi2stat,df] = chigof(observed);
        pChiDfAllResponses(pert_id,:) = [p, chi2stat,df];
    
        observed_BV = [observed(1,:);sum(observed(2:end,:))];
        [p_BV, chi2stat_BV,df_BV] = chigof(observed_BV);
        pChiDfBvResponses(pert_id,:) = [p_BV, chi2stat_BV,df_BV];
    end
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % PLOT THE RESULTS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ax1 = subplot(4,1,1);
    X = newValueMat(idx_sorted,:)';
    h1 = heatmap(X,'ColorbarVisible','off');
    h1.YDisplayLabels = param_names(pidx);
    lbls = strings(length(idx_sorted),1);
    loc_base = find(sum(X == [0;0;0;0]) == 4);
    lbls(loc_base) = "B>";
    h1.XDisplayLabels = lbls;
    title('Parameters')
    colormap(h1,[0 0.75 0.75; ...
        1 1 1; ...
        0.75 0 0.75])
    set(gca,'fontsize',10)
    
    
    subplot(4,1,[2:4])
    bar(allFrequencies(idx_sorted,:)*100,'stacked')
    colororder(sp_cols)
    xlabel('Probiotic Design')
    ylabel('Percent Samples')
    xlim([0.5 81.5])
    hold on
    xline(loc_base,'-r','LineWidth',2)
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    
    % Labels Statistics
    s = strings(size(newValueMat,1),1);
    psort = pChiDfBvResponses(idx_sorted,1)*81;
    s(psort < 0.05) = "*";
    s(psort < 0.01) = "**";
    s(psort < 0.001) = "***";
    s(psort < 0.0001) = "****";
    xticks(1:size(newValueMat,1))
    xticklabels(s)
end