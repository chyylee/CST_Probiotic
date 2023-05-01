function plot_1D_SensitivityAnalysis(all_select_outcomes,all_warnall,time_id,...
    params_to_view,param_alterkg,param_alteris,param_alterii,param_names,SSnms,sp_cols,figFlag)
    param_alter = param_alteris;
    % Loop through each parameter, plot 1D
    for param_id = 1:size(all_select_outcomes,1)
        allCounts = NaN(length(param_alter),4);
        allFrequencies = NaN(length(param_alter),4);
        for pert_id = 1:length(param_alter)
            Xwarn = squeeze(all_warnall(param_id,pert_id,:)) == 1;
            X = squeeze(all_select_outcomes(param_id,pert_id,Xwarn,time_id,2:end));
            CST_post = get_CST_post(X);
            Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
            Frequency = Counts ./ sum(Counts,2);
            allCounts(pert_id,:) = Counts;
            allFrequencies(pert_id,:) = Frequency;
        end
        
        if figFlag
            figure
        else
            subplot(2,4,param_id)
        end
        Y = allCounts;
        Y = bsxfun(@rdivide, Y, sum(Y,2))*100;
        bar(Y,'stacked')
        xticks(1:length(param_alter))
        if param_id == 1
            xticklabels(param_alterkg) % first parameter is kgrow-Probiotic (different units)
        elseif param_id == 8
            xticklabels(param_alterii)
        else
            xticklabels(param_alter)
        end 
        colororder(sp_cols)
        xlabel('Parameter Value')
        ylabel('Frequency')
        title(param_names(params_to_view(param_id)))
        xtickangle(90)
        set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
        ylim([0 100])
    end
    legend(SSnms)
end
