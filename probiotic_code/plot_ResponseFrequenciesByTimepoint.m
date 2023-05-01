function plot_ResponseFrequenciesByTimepoint(selected_timepoints,select_outcomes,warnall,SSnms,time_names,sp_cols)
    timeCounts = NaN(length(selected_timepoints),4);
    timeFrequency = NaN(length(selected_timepoints),4);
    for time_id = 1:length(selected_timepoints)
        abundance = squeeze(select_outcomes(warnall == "",selected_timepoints(time_id),2:end));
%         relative_abundance = abundance ./ sum(abundance,2);
        CST_post = get_CST_post(abundance);
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts);
        timeCounts(time_id,:) = Counts;
        timeFrequency(time_id,:) = Frequency;
    end
    
    X = timeFrequency*100;
    bar(X,'stacked')
    colororder(sp_cols)
    ylabel('Percent Samples')
    xticklabels(time_names(selected_timepoints))
    legend(SSnms)
    for i = 1:length(selected_timepoints) % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f'),'VerticalAlignment','top', 'HorizontalAlignment', 'center')
    end
    set(gca,'fontsize',10,'XColor','k','YColor','k')
end