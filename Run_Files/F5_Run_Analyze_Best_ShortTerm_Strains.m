%% F5_Run_Analyze_Best_ShortTerm_Strains.m
%
% Goal: Analyze the LHS strains that were best a promoting a certain
% response type when adminstered with a short-term probiotic regimen.
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Load Data from F4
clear;clc;
output_fdr = 'result_workspaces/F4/';
load(strcat(output_fdr,'F4_Systematic_ShortTerm_Probiotic.mat'))

result_fdr = 'result_workspaces/';
%% 2a. Run Simulation for Short Term Regimen
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
bestStrain_select_outcomes = NaN(size(bestStrains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
bestStrain_all_warnall = NaN(size(bestStrains,1),length(POPinfo.Indexes));
for combo_id = 1:size(bestStrains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    tmp = bestStrains{combo_id};
    if size(tmp,1) > 1
        sel_strain = tmp(1,:);
    else
        sel_strain = tmp;
    end
    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    bestStrain_select_outcomes(combo_id,:,:,:) = select_outcomes;
    bestStrain_all_warnall(combo_id,:) = warnall == "";
end

save(strcat(result_fdr,'F5_ShortTerm_Best_Strains.mat'),'ABXinfo','POPinfo',...
    'PROBinfo','bestStrain_all_warnall','bestStrain_select_outcomes','-append')

%% 2b. Plot Result for Short Term Regimen

% pulls the probiotic parameters
definedStrain = NaN(size(bestStrains,1),8);
for k = 1:size(bestStrains,1)
    tmp = bestStrains{k};
    if size(tmp,1) > 1
        sel_strain = tmp(1,:);
    else
        sel_strain = tmp;
    end
    definedStrain(k,:) = sel_strain;
end

doseLabels = PROBinfo.customProbDosing-1;
plotLabels = {'Promotes nAB dominance', 'Promotes Li dominance',...
    'Promotes oLB dominance','Promotes P dominance'};

plot_BestLHS_Strains(PROBinfo,ABXinfo,bestStrain_select_outcomes,...
    bestStrain_all_warnall,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 8])


%% 3a. Run Simulation with Lactin-V Regimen
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Lactin-V Regimen
simulationType = "comboLV";
ws_name = 'input_files/SSConfig-Analysis-HMP-Virtual-Population.mat';
[POPinfo,PROBinfo,ABXinfo,~] = define_Common_Simulation_Inputs(ws_name,simulationType);
ABXinfo.dose = -2.64;

PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
bestStrain_select_outcomes_LV = NaN(size(bestStrains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
bestStrain_all_warnall_LV = NaN(size(bestStrains,1),length(POPinfo.Indexes));
for combo_id = 1:size(bestStrains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    tmp = bestStrains{combo_id};
    if size(tmp,1) > 1
        sel_strain = tmp(1,:);
    else
        sel_strain = tmp;
    end
    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    bestStrain_select_outcomes_LV(combo_id,:,:,:) = select_outcomes;
    bestStrain_all_warnall_LV(combo_id,:) = warnall == "";
end

save(strcat(result_fdr,'F5_ShortTerm_Best_Strains.mat'),...
    'bestStrain_all_warnall_LV','bestStrain_select_outcomes_LV','-append')

%% 3b. Plot Resuls for Lactin-V Regimen

% pulls the probiotic parameters
definedStrain = NaN(size(bestStrains,1),8);
for k = 1:size(bestStrains,1)
    tmp = bestStrains{k};
    if size(tmp,1) > 1
        sel_strain = tmp(1,:);
    else
        sel_strain = tmp;
    end
    definedStrain(k,:) = sel_strain;
end

lvPROBinfo = PROBinfo;
lvABXinfo = ABXinfo;
lvABXinfo.addABX = true; % true: administer ABX, false: no ABX

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([0 1 0 0 1 0 0],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values
lvPROBinfo.customProbDosing = customDosing;

doseLabels = lvPROBinfo.customProbDosing+4;
plotLabels = {'Promotes nAB dominance', 'Promotes Li dominance',...
    'Promotes oLB dominance','Promotes P dominance'};

plot_BestLHS_Strains(lvPROBinfo,lvABXinfo,bestStrain_select_outcomes_LV,...
    bestStrain_all_warnall_LV,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 8])

