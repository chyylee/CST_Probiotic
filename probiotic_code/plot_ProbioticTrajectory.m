function plot_ProbioticTrajectory(select_plot,eval_points,sp_cols)
    spl_sz = size(select_plot,1);
    comb_avg = [squeeze(mean(select_plot,1))];
    comb_std = [squeeze(std(select_plot,[],1))];
    comb_CI = tinv(0.99,spl_sz-1)*comb_std/sqrt(spl_sz);
    upper_CI = comb_avg + comb_CI;
    lower_CI = comb_avg - comb_CI;

    for i = 1:4
        curve1 = upper_CI(:,i)';
        curve2 = lower_CI(:,i)';
        inBetweenRegionX = [eval_points, fliplr(eval_points)];
        inBetweenRegionY = [curve1, fliplr(curve2)];
        fill(inBetweenRegionX, inBetweenRegionY, sp_cols(i,:),...
            'FaceAlpha',0.3,'EdgeColor',sp_cols(i,:));
        hold on
        p = plot(eval_points,comb_avg(:,i),'LineWidth',1.5,'Color',sp_cols(i,:));
    end
end