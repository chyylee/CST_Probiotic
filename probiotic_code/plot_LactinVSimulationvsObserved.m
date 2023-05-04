function plot_LactinVSimulationvsObserved(fdr_nm,ws_list)
    [~,~,sp_cols,~] = get_naming_terms();
    % FORMAT DATA & STATISTICS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    % ~~~~~~~~~~~~ Cohen et al 2020 Results ~~~~~~~~~~~ 
    ClinicalRecur = [46 59; % count @ 12 wk, 24 wk Lactin-V
        34 41]; % count@ 12 wk, 24 wk Lactin-V
    ClinicalCure = [87 63;
        30 21];
    ClinicalN = ClinicalCure + ClinicalRecur;
    ClinicalResults = (ClinicalRecur./ClinicalN)*100;
    
    % ~~~~~~~~~~~~ Pull Model Results ~~~~~~~~~~~ 
    null = load(strcat(fdr_nm,'F3_','Null_LactinVPlacebo.mat'));
    
    % ~~~~~~~~~~~~ COMPARE PLACEBO ARM TO STRAINS 12 WKS ~~~~~~~~~~~ 
    abundance = squeeze(null.all_select_outcomes(2,null.all_warnall(2,:) == 1,12,2:end));
    CST_post = get_CST_post(abundance);
    Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
    Frequency = Counts ./ sum(Counts) *100;
    placeboClin = [ClinicalResults(2,1);
        100 - ClinicalResults(2,1);
        0;
        0];
    
    % ~~~~~~~~~~~~ 12WKS: COMPARE PLACEBO ARM TO STRAINS 12 WKS ~~~~~~~~~~~ 
    subplot(2,3,1)
    n = 2;
    X = [placeboClin';Frequency];
    N = ClinicalN(2,1);
    Y = [X(1,1)/100*N, N - X(1,1)/100*N;
        X(2,1)/100*sum(Counts), sum(Counts) - X(2,1)/100*sum(Counts)];
    [pPLCBO12,~,~] = chigof(Y');
    superbar(X(:,1),'P',[NaN,pPLCBO12;pPLCBO12,NaN],'BarFaceColor', sp_cols(1,:),...
        'BarEdgeColor', 'k','PStarFontSize',10,'BarLineWidth',1); % placebo 12 wk
    for i = 1:n % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f%%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center')
    end
    text(1.5,max(X(:,1))+15,strcat("P = ",num2str(pPLCBO12,'%.4f')),"HorizontalAlignment","center")
    xticklabels({'Clinical','Model'}); ylabel('Percent Samples'); 
    title('12 weeks, Placebo'); ylim([0 100]); xticks(1:2);
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    
    
    % ~~~~~~~~~~~~ 24WKS: COMPARE PLACEBO ARM TO STRAINS 24 WKS ~~~~~~~~~~~ 
    abundance = squeeze(null.all_select_outcomes(2,null.all_warnall(2,:) == 1,14,2:end));
    CST_post = get_CST_post(abundance);
    Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
    Frequency = Counts ./ sum(Counts) *100;
    
    placeboClin = [ClinicalResults(2,2);
        100 - ClinicalResults(2,2);
        0;
        0];
    
    subplot(2,3,4)
    X = [placeboClin';Frequency];
    N = ClinicalN(2,2);
    Y = [X(1,1)/100*N, N - X(1,1)/100*N;
        X(2,1)/100*sum(Counts), sum(Counts) - X(2,1)/100*sum(Counts)];
    [pPLCBO12,~,~] = chigof(Y');
    superbar(X(:,1),'P',[NaN,pPLCBO12;pPLCBO12,NaN],'BarFaceColor', sp_cols(1,:),...
        'BarEdgeColor', 'k','PStarFontSize',10,'BarLineWidth',1); % placebo 12 wk
    for i = 1:n % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f%%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center')
    end
    title('24 weeks, Placebo'); ylim([0 100]); xticklabels({'Clinical','Model'})
    ylabel('Percent Samples'); xticks(1:2);
    text(1.5,max(X(:,1))+15,strcat("P = ",num2str(pPLCBO12,'%.4f')),"HorizontalAlignment","center")
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % LACTIN-V ARM
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LVClin12 = [ClinicalResults(1,1);
        100 - ClinicalResults(1,1);
        0;
        0];
    
    LVClin24 = [ClinicalResults(1,2);
        100 - ClinicalResults(1,2);
        0;
        0];
    
    N12 = ClinicalN(1,1);
    N24 = ClinicalN(1,2);
    LVClin = [LVClin12(1,1)/100*N12, N12 - LVClin12(1,1)/100*N12;
        LVClin24(1,1)/100*N24, N24 - LVClin24(1,1)/100*N24];
    
    ws_list_loc = strcat(fdr_nm,ws_list);
    all_Counts = NaN(2,length(ws_list),4);
    PvsC = NaN(2,length(ws_list));
    wk_indx = [12 14];
    for wk = 1:length(wk_indx)
        for i = 1:length(ws_list)
            ws = load(ws_list_loc{i});
            abundance = squeeze(ws.all_select_outcomes(1,ws.all_warnall(1,:) == 1,wk_indx(wk),2:end));
            CST_post = get_CST_post(abundance);
            Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
            all_Counts(wk,i,:) = Counts;
            Y = [LVClin(wk,:); 
                Counts(1), sum(Counts) - Counts(1)];
            [PvsC(wk,i),~,~] = chigof(Y');
        end
    end
    Frequency = squeeze(all_Counts(1,:,:)) ./ sum(squeeze(all_Counts(1,:,:)),2)*100;
    
    % ~~~~~~~~~~~~ 12WKS: COMPARE LV ARM TO STRAINS 12 WKS ~~~~~~~~~~~ 
    subplot(2,3,[2,3])
    n = 1+size(Frequency,2);
    X = [LVClin12';Frequency];
    P = NaN(n,n);
    P(1,2:end) = PvsC(1,:);
    P(2:end,1) = PvsC(1,:);
    superbar(X(:,1),'P',P,'BarFaceColor', sp_cols(1,:),...
        'BarEdgeColor', 'k','PStarFontSize',10,'BarLineWidth',1); % placebo 12 wk
    for i = 1:n % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f%%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center')
    end
    str = 'P_1 = %.4f; P_2 = %.4f; P_3 = %.4f; P_4 = %.4f';
    text(3,110,sprintf(str,PvsC(1,:)),"HorizontalAlignment","center"); xticks(1:n);
    ylim([0 120]); xticklabels(horzcat({'Clinical'},extractBetween(ws_list,'F3_','_')))
    title('12 weeks, Treatment'); ylabel('Percent Samples')
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    
    % ~~~~~~~~~~~~ 24WKS: COMPARE LV ARM TO STRAINS 24 WKS ~~~~~~~~~~~ 
    Frequency = squeeze(all_Counts(2,:,:)) ./ sum(squeeze(all_Counts(2,:,:)),2)*100;
    subplot(2,3,[5,6])
    n = 1+size(Frequency,2);
    X = [LVClin24';Frequency];
    P = NaN(n,n);
    P(1,2:end) = PvsC(2,:);
    P(2:end,1) = PvsC(2,:);
    superbar(X(:,1),'P',P,'BarFaceColor', sp_cols(1,:),...
        'BarEdgeColor',  'k','PStarFontSize',10,'BarLineWidth',1); % placebo 12 wk
    for i = 1:n % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f%%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center')
    end
    str = 'P_1 = %.4f; P_2 = %.4f; P_3 = %.4f; P_4 = %.4f';
    text(3,110,sprintf(str,PvsC(2,:)),"HorizontalAlignment","center"); xticks(1:n);
    ylim([0 120]); xticklabels(horzcat({'Clinical'},extractBetween(ws_list,'F3_','_L')))
    title('24 weeks, Treatment'); ylabel('Percent Samples')
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
end