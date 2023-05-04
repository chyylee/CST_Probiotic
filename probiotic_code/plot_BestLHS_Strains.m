function plot_BestLHS_Strains(PROBinfo,ABXinfo,bestStrain_select_outcomes, bestStrain_all_warnall,definedStrain,doseLabels,plotLabels)
    [param_names,SSnms,sp_cols,~] = get_naming_terms();
    
    numTimepoints = length(PROBinfo.customEvalPoints);
    numStrains = size(bestStrain_select_outcomes,1);
    
    
    for id = 1:numStrains
        select_outcomes = squeeze(bestStrain_select_outcomes(id,bestStrain_all_warnall(id,:)==1,1:numTimepoints ,:));
        select_plot= NaN(size(select_outcomes,1),numTimepoints,4);
        for i = 1:size(select_outcomes)
            select_plot(i,:,:) = squeeze(select_outcomes(i,:,2:end)) ./ sum(squeeze(select_outcomes(i,:,2:end)),2);
        end
        eval_points = PROBinfo.customEvalPoints;
        
        
        % ~~~~~~ Network
        subplot(numStrains,12,[1 2] + (id-1)*12)
        sp_names = {'nAB','Li','oLB','Probiotic'};
        selectStrain = definedStrain(id,:);
        ttl = '';
        pidx = PROBinfo.strainCharIndexes;
        line_scale = 1e-20; % line thickness 
        % format probiotic strain paarameters for plotting
        tmp_net_vect = line_scale*ones(size(param_names));
        tmp_net_vect(pidx) = selectStrain;
        tmp_net_vect(tmp_net_vect == 0) = line_scale;
        % Get Thicknesses
        tmp_thickness_vect = line_scale*ones(size(param_names));
        tmp_thickness_vect(pidx) = selectStrain;
        tmp_thickness_vect(tmp_thickness_vect == 0) = line_scale;
        sclThck = 3;
        net_vect = tmp_net_vect(5:end)';
        intDir = net_vect;
        intThick = abs(sclThck*tmp_thickness_vect(5:end));
        plot_network(net_vect,sp_names,sp_cols,intDir,intThick,ttl)
        caxis([-0.12,0.12])
        colorbar
        title(strcat(plotLabels{id}, " Characteristics"))
        
        % ~~~~~~ Abundance-Time Profile
        subplot(numStrains,12,[4 5 6 7] + (id-1)*12)
        plot_ProbioticTrajectory(select_plot,eval_points+5,sp_cols)
        xlim([-5,60])
        ylim([0 1]); ylabel('Relative Abundance');xlabel('time (d)')
        xline(doseLabels,':',{'Dose'})
    
        if ABXinfo.addABX
            patch([0 5 5 0], [0 0 1 1], [1 0 0], 'FaceAlpha', 0.1);
        end
    
        set(gca,'fontsize',10,'XColor','k','YColor','k')
        title(strcat("Average Abundance-Time Profile ", plotLabels{id}))
        
        warnall = strings(1,2000);
        warnall(squeeze(bestStrain_all_warnall(id,:,:,:))==0) = "wng";
        
        select_outcomes = squeeze(bestStrain_select_outcomes(id,:,:,:));
        selected_timepoints = [1,numTimepoints];
        time_names = string(eval_points);
        
        % ~~~~~~ BV Recurrence by 12mo
        subplot(numStrains,12,[9 10 11]+ (id-1)*12)
        [Counts] = plot_ResponseFrequenciesByTimepoint(selected_timepoints,select_outcomes,warnall,SSnms,time_names,sp_cols);
        xlim([1.5 2.5]); xticks(1:2);xticklabels({'','12mo'});
        hold on
        X = Counts ./ sum(Counts,2)*100;
        Y = cumsum(X,2);
        i = 2;
        
        handlesToChildObjectsInAxes = findobj(gca, 'Type', 'text');
        if ~isempty(handlesToChildObjectsInAxes)
            delete(handlesToChildObjectsInAxes);
        end
        
        for j = 1:4
            if j ~= 4 && X(i,j) > 0
                text(i,Y(i,j),num2str(X(i,j),'%.1f %%'),'VerticalAlignment','top', 'HorizontalAlignment',...
                    'center','Color','k')
            elseif j == 4 && X(i,j) > 0
                text(i,Y(i,j),num2str(X(i,j),'%.1f %%'),'VerticalAlignment','top', 'HorizontalAlignment',...
                    'center','Color','w')
            else
                text(i,0,"")
            end
        
        end
        legend('Location','eastoutside')
        title("12mo status")
    end
end