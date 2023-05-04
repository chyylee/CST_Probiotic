%% F5_Run_Analyze_LV_Dosing_Alernatives.m
%
% Goal: Visualize how different dosing strategies impact Lactin-V efficacy
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Load Data
load('F5_LV_Best_Strains.mat')

%% 2. Pull 12 week and 24 week data

%  ~~~~~~~~~~~~ Cohen et al 2020 Results ~~~~~~~~~~~ 
ClinicalRecur = [46 59; % count @ 12 wk, 24 wk Lactin-V
    34 41]; % count@ 12 wk, 24 wk Lactin-V
ClinicalCure = [87 63;
    30 21];
ClinicalN = ClinicalCure + ClinicalRecur;
ClinicalResults = (ClinicalRecur./ClinicalN)*100;
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   

% select_outcomes: 5 x 2000 x 89 x 5
    % Index 78 = 12 wk, Index 80 = 24 wk

LV_Output = {LVStrain_select_outcomes_biwk,... % 1 dose / 2 weeks
    LVStrain_select_outcomes_1wk,... % 1 dose / 1 week
    LVStrain_select_outcomes_2wk,... % 2 doses / 1 week
    LVStrain_select_outcomes_4wk,... % 4 doses / 1 week
    LVStrain_select_outcomes_7wk}; % daily

LV_Warning = {LVStrain_all_warnall_biwk,... % 1 dose / 2 weeks
    LVStrain_all_warnall_1wk,... % 1 dose / 1 week
    LVStrain_all_warnall_2wk,... % 2 doses / 1 week
    LVStrain_all_warnall_4wk,... % 4 doses / 1 week
    LVStrain_all_warnall_7wk}; % daily

wk_indx = [78 80]; % 12 wk, 24 wk

all_Counts = NaN(5,length(wk_indx),length(LV_Output),4);
all_Frequency = NaN(5,length(wk_indx),length(LV_Output),4);

for wk = 1:length(wk_indx)
    for reg_id = 1:length(LV_Output)
        atmp_all_select_outcomes = LV_Output{reg_id};
        atmp_all_warnall = LV_Warning{reg_id};
        for strain_id = 1:5
            tmp_all_select_outcomes = squeeze(atmp_all_select_outcomes(strain_id,:,:,:));
            tmp_all_warnall = atmp_all_warnall(strain_id,:);
            abundance = squeeze(tmp_all_select_outcomes(tmp_all_warnall == 1,wk_indx(wk),2:end));
            CST_post = get_CST_post(abundance);
            Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
            all_Counts(strain_id,wk,reg_id,:) = Counts;
            all_Frequency(strain_id,wk,reg_id,:) = Counts ./ sum(Counts,2);
        end
    end
end

% Counts = 5 x 2 x 5 x 4 (strain x week x regimen x response)

%% 3. Plot by week & Strain
resp_id = 1; % nAB dominant
wk = 1;

subplot(2,1,1)
wk12_data = squeeze(all_Frequency(:,wk,:,resp_id))*100;
b = bar(wk12_data);
strainLabels = "LV Strain #" + string(1:5);
set(gca,'xtick',[1:5],'xticklabels',strainLabels)
colororder(parula(5))
title('Week 12'); ylabel('Failure Rate')
p = yline(ClinicalResults(1,1),":",{'Clinical Observation'});
legend(b,{'Bi-weekly','weekly', '2x weekly','4x weekly','Daily'},'location','southoutside',...
    'Orientation','horizontal')

subplot(2,1,2)
wk = 2;
wk12_data = squeeze(all_Frequency(:,wk,:,resp_id))*100;
b = bar(wk12_data);
strainLabels = "LV Strain #" + string(1:5);
set(gca,'xtick',[1:5],'xticklabels',strainLabels)
colororder(parula(5))
title('Week 24'); ylabel('Failure Rate')
p = yline(ClinicalResults(1,2),":",{'Clinical Observation'});
legend(b,{'Bi-weekly','weekly', '2x weekly','4x weekly','Daily'},'location','southoutside',...
    'Orientation','horizontal')

%% 4. Combine strains, compare regimens
resp_id = 1; % nAB dominant
ttls = {'Week 12','Week 24'};

c = [1 1 1;
    1 1 1;
    0 0 0;
    0 0 0;
    0 0 0;
    0 0 0];

for wk = 1:2
    numRegimens = size(all_Counts,3);
    for i = 1:numRegimens
        for j = 1:numRegimens
            tmp1 = squeeze(all_Frequency(:,wk,i,resp_id))*100;
            tmp2 = squeeze(all_Frequency(:,wk,j,resp_id))*100;
            [p,h] = ranksum(tmp1,tmp2);
            P(i,j) = p;
        end
    end
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P);

    average = mean(squeeze(all_Frequency(:,wk,:,resp_id)*100));
    std_dev = std(squeeze(all_Frequency(:,wk,:,resp_id)*100));
    
    subplot(1,2,wk)
    superbar(average,'E',std_dev,'P',P*numRegimens,'BarFaceColor',...
        brighten(parula(numRegimens ),-0.25),'BarEdgeColor','k',...
        'BarLineWidth',1,'PStarOffset',2.25,'PStarFontSize',8,...
        'PLineWidth',0.75)
    set(gca,'xtick',[1:numRegimens],'xticklabels',{'Bi-weekly',...
        'Weekly', '2x weekly','4x weekly','Daily'})
    hold on
    swarmchart(repmat(1:5,5,1),squeeze(all_Frequency(:,wk,:,resp_id)*100),...
        40, brighten(parula(numRegimens),-0.9),'filled','XJitterWidth',0.5)
    ylabel('Failure Rate'); title(ttls{wk}); ylim([0 100])
    set(gca,'fontsize',8)
    hold on
    for i = 1:length(average) % Label for percent nAB response post
        text(i,10,num2str(average(i),'%.1f %%'),'VerticalAlignment','top',...
            'HorizontalAlignment', 'center','FontSize',8,'Color',c(i,:))
    end
end
set(gcf,'units','inches','position',[0 0 12 3])

