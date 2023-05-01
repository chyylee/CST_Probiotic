%% plot_trajectories_indv(SSnms, CST_post, warnall, tall, yall, sp_cols, tspanP)
%
% Use: Visualize trajectories for a certain CST_post ABX
%
% SSnms: cell array of equilibrium state names
% CST_post: numerical array of CST post probiotic (index matches order in
%           SSnms)
% warnall: probiotic runs that had a ODE warning thrown
% tall: cell array of the time steps of the ode sim
% yall: cell array of the absolute abundance of each microbe in the ode sim
%       over time
% sp_cols: colors for plot
% tspanP: when the probiotic was added
%
%%
function plot_trajectories_indv(SSnms, CST_post, warnall, tall, yall, sp_cols, dose_lb)
    [indx,~] = listdlg('ListString',SSnms);

    sidx = find(CST_post' == indx & warnall == "");
    
    if length(sidx) > 40
        L = 40;
    else
        L = length(sidx);
    end

    m = ceil(sqrt(L));
    n = ceil(L/m);
    for i = 1:L
        ttmp = tall{sidx(i)};
        ytmp = yall{sidx(i)};

        figure(i)
%         subplot(m,n,i)
        plot(ttmp,ytmp./sum(ytmp,2),'linewidth',3)
        hold on
        for k = 1:length(dose_lb)
            xline(dose_lb(k),':')
        end
        xlabel('time')
        ylabel('Relative Abundance')
        title(strcat("Set #: ", num2str(sidx(i))))
        set(gca,'fontsize',14)
        xlim([10,108])
        ylim([0,1])

        colororder(sp_cols)
    end
    legend({'NO','Li','oLB','Probiotic','Dose Point'})
    
end