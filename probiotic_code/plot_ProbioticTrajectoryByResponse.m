function plot_ProbioticTrajectoryByResponse(resp_id,select_outcomes,warnall,tall,yall,time_id,dose_lbls,SSnms,sp_cols)
    abundance = squeeze(select_outcomes(:,time_id,2:end));
    CST_post = get_CST_post(abundance);
    c = 1;

    resp_indx = CST_post == resp_id & warnall' == "";
    tsub = tall(resp_indx);
    ysub = yall(resp_indx);
    % Plot the average and std
    eval_points = [0:0.5:30,31:10:380]; % pull subset of time points
    select_plot = NaN(length(tsub),length(eval_points),4);
    for net_id = 1:length(tsub)
       tcol = tsub{net_id};
       ycol = ysub{net_id};
       for i = 1:length(eval_points)
            [~,idx] = min(abs(tcol - eval_points(i))); % finds closest time point
            if ~isempty(idx)
                select_plot(net_id,i,:) = ycol(idx(1),:) ./ sum(ycol(idx(1),:),2);
            end
       end
    end
      
    plot_ProbioticTrajectory(select_plot,eval_points,sp_cols)
    xline(dose_lbls,':',{'Dose'}); xline(dose_lbls(end)+30,'-',{'1 mo'})
    xlim([0 120]); xlabel('Days'); ylabel('Relative Abundance')
    title(strcat("Response of: "), strcat(SSnms(resp_id), " N = ", num2str(sum(resp_indx))))
    set(gca,'fontsize',10,'XColor','k','YColor','k')
    c = c + 1;

end