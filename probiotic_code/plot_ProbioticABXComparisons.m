function [Counts,all_PvaluesLB,all_PvaluesBV] = ...
    plot_ProbioticABXComparisons(all_select_outcomes,all_warnall,...
    testDesigns,testDesignNames,timePoints)

    % FORMAT DATA & STATISTICS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    all_PvaluesLB = NaN(length(timePoints),length(testDesignNames),length(testDesignNames));
    all_PvaluesBV = NaN(length(timePoints),length(testDesignNames),length(testDesignNames));
    Frequencies = NaN(length(timePoints),length(testDesigns),4);
    Counts = NaN(length(timePoints),length(testDesigns),4);
    c = 1;
    for time_id = timePoints    
        for test_id = 1:length(testDesignNames)
            nonErrRuns = all_warnall(test_id,:) == 1;
            selOutcome = squeeze(all_select_outcomes(test_id,nonErrRuns,time_id,2:end));
            CST_post = get_CST_post(selOutcome);
            tmp_Counts = [sum(CST_post == 1);sum(CST_post == 2);sum(CST_post == 3);sum(CST_post == 4)];
            Frequencies(c,test_id,:) = tmp_Counts ./ sum(tmp_Counts);
            Counts(c,test_id,:) = tmp_Counts;
        end
    
        C = squeeze(Counts(c,:,:));
        P = NaN(size(C,1),size(C,1));
        for i = 1:size(C,1)
            for j = 1:size(C,1)
                observed = [C(i,:)', C(j,:)'];
                % remove entries with zeros
                tmp = sum(observed == 0,2);
                if sum(tmp) == 2
                    observed = observed(1:end-1,:);
                end
                [p, chi2stat,df] = chigof(observed);
                P(i,j) = p;
                observed = [C(i,:)', C(j,:)'];
                observedLB = [observed(1,:); sum(observed(2:end,:))];
                [p, chi2stat,df] = chigof(observedLB);
                PLB(i,j) = p;
            end
        end
        all_PvaluesBV(c,:,:) = P;   
        all_PvaluesLB(c,:,:) = PLB;  
        c = c + 1;
    end
    
    
    
    % PLOT THE RESULTS 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    [~,~,sp_cols,time_names] = get_naming_terms;
    sel_time = timePoints;
    for time_id = [1 2 3 4]
        subplot(1,4,time_id)
        Freq = squeeze(Frequencies(time_id,:,:))*100;
        P = round(squeeze(all_PvaluesLB(time_id,:,:)),6);
        bar(Freq,'stacked')
        set(gca,'Xticklabels',{'Pr','Pr+ABX','ABX'})
        xtickangle(90)
        colororder(sp_cols)
        hold on
    
        Y = ones(length(testDesigns),1)*100;
        superbar(Y,'P',P,'BarFaceColor', 'none', 'BarEdgeColor', 'none',...
            'PStarOffset',4,'PLineWidth',0.75)
        for i = 1:size(Y,1) % Label for percent nAB response post
            text(i,Freq(i,1),num2str(Freq(i,1),'%.1f%%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center',...
                'fontname','Arial','Color','k')
        end
        title(time_names(sel_time(time_id)))
        ylabel('Percent Samples')
        set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
        ylim([0 140])
    end
end

