%% F5_Run_Analyze_Best_LactinV_Strains.m
%
% Goals:
%   1) Visualize strains that promote nAB-, Li-, oLB-, P-dominance
%       consistently across the virtual population
%   2) Test probiotic regimens for possible Lactin-V strains
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. Load Data (From F4)
clear;clc;
output_fdr = 'result_workspaces/F4/';
load(strcat(output_fdr,'F4_Systematic_LactinV_Probiotic.mat'))
%   - bestStrains (cell array with the strain that (1) promoted the highest
%   frequency of nAB dominated responses across patients (2) Li dominated
%   (3) oLB dominated (4) Probiotic dominated (5) candidate Lactin-V
%   strains

result_fdr = 'result_workspaces/';

%% 2a. Run Simulation LV Simulation with Best Strains

PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
bestStrain_select_outcomes = NaN(size(bestStrains,1),length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
bestStrain_all_warnall = NaN(size(bestStrains,1),length(POPinfo.Indexes));
definedStrain = NaN(size(bestStrains,1),8);
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
    %definedStrain(combo_id,:) = sel_train;
end

save(strcat(result_fdr,'F5_LV_Best_Strains.mat'),'ABXinfo','POPinfo',...
    'PROBinfo','bestStrain_all_warnall','bestStrain_select_outcomes')

%% 2b. Plot Result best Strains

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

doseLabels = PROBinfo.customProbDosing+4;
plotLabels = {'Promotes nAB dominance', 'Promotes Li dominance',...
    'Promotes oLB dominance','Promotes P dominance', 'Best Fit LV Strain'};

plot_BestLHS_Strains(PROBinfo,ABXinfo,bestStrain_select_outcomes,...
    bestStrain_all_warnall,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])


%% 3a. Evaluate LV strains with Short-term dosing
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wk = [1 1 1 1 1 1 1]; % wk 1: 1/day 4 days
customDosing = find(wk); % get numerical values
PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;

ABXinfo.addABX = false; % true: administer ABX, false: no ABX

PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
bestStrain_select_outcomes_ST = NaN(size(bestStrains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
bestStrain_all_warnall_ST = NaN(size(bestStrains,1),length(POPinfo.Indexes));
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

    bestStrain_select_outcomes_ST(combo_id,:,:,:) = select_outcomes;
    bestStrain_all_warnall_ST(combo_id,:) = warnall == "";
end

save(strcat(output_fdr_nm,'F5_LV_Best_Strains.mat'),...
    'bestStrain_all_warnall_ST','bestStrain_select_outcomes_ST','-append')

%% 3a. Plot Result Short-term dosing

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

% re-define simulations entries
stPROBinfo = PROBinfo;
wk = [1 1 1 1 1 1 1]; % wk 1: 1/day 4 days
customDosing = find(wk); % get numerical values
stPROBinfo.custom = true;
stPROBinfo.customProbDosing = customDosing;
stABXinfo.addABX = false;

doseLabels = PROBinfo.customProbDosing-1;
plotLabels = {'Promotes nAB dominance', 'Promotes Li dominance',...
    'Promotes oLB dominance','Promotes P dominance', 'Best Fit LV Strain'};

plot_BestLHS_Strains(stPROBinfo,stABXinfo,bestStrain_select_outcomes_ST,...
    bestStrain_all_warnall_ST,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])


%% 4a. LV Strains as Acute Dosing
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wk = [1 1 1 1 1 1 1]; % wk 1: 1/day 4 days
customDosing = find(wk); % get numerical values
PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;

ABXinfo.addABX = false; % true: administer ABX, false: no ABX

PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
LVStrain_select_outcomes = NaN(size(LV_strains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
LVStrain_all_warnall = NaN(size(LV_strains,1),length(POPinfo.Indexes));
for combo_id = 1:size(LV_strains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    sel_strain = LV_strains(combo_id,:);

    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    LVStrain_select_outcomes(combo_id,:,:,:) = select_outcomes;
    LVStrain_all_warnall(combo_id,:) = warnall == "";
end

save(strcat(output_fdr_nm,'F5_LV_Best_Strains.mat'),...
    'LVStrain_all_warnall','LVStrain_select_outcomes','-append')


%% 4b. Plot Results LV Strains Acute Dosing

stPROBinfo = PROBinfo;
wk = [1 1 1 1 1 1 1]; % wk 1: 1/day 4 days
customDosing = find(wk); % get numerical values
stPROBinfo.custom = true;
stPROBinfo.customProbDosing = customDosing;
stABXinfo.addABX = false;

doseLabels = PROBinfo.customProbDosing-1;
definedStrain = LV_strains;

plotLabels = cellstr("LV Strain #" + string(1:5));
plot_BestLHS_Strains(stPROBinfo,stABXinfo,LVStrain_select_outcomes,...
    LVStrain_all_warnall,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])

%% 5a. LV Strains as LV Regimen (1/wk)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ABXinfo.addABX = true; % true: administer ABX, false: no ABX
ABXinfo.dose = -2.64; % decay rate
ABXinfo.duration = 5; % days

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([0 0 0 1 0 0 0],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;
PROBinfo.numberDoses = 1; % days
PROBinfo.doseFrequency = NaN; % every XX days


PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
LVStrain_select_outcomes_1wk = NaN(size(LV_strains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
LVStrain_all_warnall_1wk = NaN(size(LV_strains,1),length(POPinfo.Indexes));
for combo_id = 1:size(LV_strains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    sel_strain = LV_strains(combo_id,:);

    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    LVStrain_select_outcomes_1wk(combo_id,:,:,:) = select_outcomes;
    LVStrain_all_warnall_1wk(combo_id,:) = warnall == "";
end

save(strcat(result_fdr,'/F5/F5_LV_Best_Strains.mat'),...
    'LVStrain_all_warnall_1wk','LVStrain_select_outcomes_1wk','-append')


%% 5b. Plot Results LV LV Regimen (1/wk)

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([0 0 0 1 0 0 0],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

lvPROBinfo = PROBinfo;
lvPROBinfo.customProbDosing = customDosing;

doseLabels = lvPROBinfo.customProbDosing+4;
definedStrain = LV_strains;

plotLabels = cellstr("LV Strain #" + string(1:5));
plot_BestLHS_Strains(lvPROBinfo,ABXinfo,LVStrain_select_outcomes_1wk,...
    LVStrain_all_warnall_1wk,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])

%% 6a. LV Strains as LV Regimen (2/wk)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ABXinfo.addABX = true; % true: administer ABX, false: no ABX
ABXinfo.dose = -2.64; % decay rate
ABXinfo.duration = 5; % days

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([0 0 1 0 1 0 0],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;
PROBinfo.numberDoses = 1; % days
PROBinfo.doseFrequency = NaN; % every XX days


PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
LVStrain_select_outcomes_2wk = NaN(size(LV_strains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
LVStrain_all_warnall_2wk = NaN(size(LV_strains,1),length(POPinfo.Indexes));
for combo_id = 1:size(LV_strains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    sel_strain = LV_strains(combo_id,:);

    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    LVStrain_select_outcomes_2wk(combo_id,:,:,:) = select_outcomes;
    LVStrain_all_warnall_2wk(combo_id,:) = warnall == "";
end

save(strcat(result_fdr,'/F5/F5_LV_Best_Strains.mat'),...
    'LVStrain_all_warnall_2wk','LVStrain_select_outcomes_2wk','-append')


%% 6b. Plot Results LV Strains LV Regimen (2/wk)

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([0 0 1 0 1 0 0],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

lvPROBinfo = PROBinfo;
lvPROBinfo.customProbDosing = customDosing;

doseLabels = lvPROBinfo.customProbDosing+4;
definedStrain = LV_strains;

plotLabels = cellstr("LV Strain #" + string(1:5));
plot_BestLHS_Strains(lvPROBinfo,ABXinfo,LVStrain_select_outcomes_2wk,...
    LVStrain_all_warnall_2wk,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])


%% 7a. LV Strains Bi-weekly
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ABXinfo.addABX = true; % true: administer ABX, false: no ABX
ABXinfo.dose = -2.64; % decay rate
ABXinfo.duration = 5; % days

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([0 0 0 0 0 0 1 0 0 0 0 0 0],1,wks/2);
customDosing = find([wk1,int_lbs]); % get numerical values

PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;
PROBinfo.numberDoses = 1; % days
PROBinfo.doseFrequency = NaN; % every XX days


PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
LVStrain_select_outcomes_biwk = NaN(size(LV_strains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
LVStrain_all_warnall_biwk = NaN(size(LV_strains,1),length(POPinfo.Indexes));
for combo_id = 1:size(LV_strains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    sel_strain = LV_strains(combo_id,:);

    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    LVStrain_select_outcomes_biwk(combo_id,:,:,:) = select_outcomes;
    LVStrain_all_warnall_biwk(combo_id,:) = warnall == "";
end

save(strcat(result_fdr,'/F5/F5_LV_Best_Strains.mat'),...
    'LVStrain_all_warnall_biwk','LVStrain_select_outcomes_biwk','-append')


%% 7b. Plot Results LV Strains Bi-weekly

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([0 0 0 0 0 0 1 0 0 0 0 0 0],1,wks/2);
customDosing = find([wk1,int_lbs]); % get numerical values

lvPROBinfo = PROBinfo;
lvPROBinfo.customProbDosing = customDosing;

doseLabels = lvPROBinfo.customProbDosing+4;
definedStrain = LV_strains;

plotLabels = cellstr("LV Strain #" + string(1:5));
plot_BestLHS_Strains(lvPROBinfo,ABXinfo,LVStrain_select_outcomes_biwk,...
    LVStrain_all_warnall_biwk,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])

%% 8a. LV Strains every 4x week
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ABXinfo.addABX = true; % true: administer ABX, false: no ABX
ABXinfo.dose = -2.64; % decay rate
ABXinfo.duration = 5; % days

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([1 0 1 0 1 0 1],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;
PROBinfo.numberDoses = 1; % days
PROBinfo.doseFrequency = NaN; % every XX days


PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
LVStrain_select_outcomes_4wk = NaN(size(LV_strains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
LVStrain_all_warnall_4wk = NaN(size(LV_strains,1),length(POPinfo.Indexes));
for combo_id = 1:size(LV_strains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    sel_strain = LV_strains(combo_id,:);

    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    LVStrain_select_outcomes_4wk(combo_id,:,:,:) = select_outcomes;
    LVStrain_all_warnall_4wk(combo_id,:) = warnall == "";
end

save(strcat(result_fdr,'/F5/F5_LV_Best_Strains.mat'),...
    'LVStrain_all_warnall_4wk','LVStrain_select_outcomes_4wk','-append')


%% 8b. Plot Results LV 4x

% probiotic dosing regimen
wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([1 0 1 0 1 0 1],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

lvPROBinfo = PROBinfo;
lvPROBinfo.customProbDosing = customDosing;

doseLabels = lvPROBinfo.customProbDosing+4;
definedStrain = LV_strains;

plotLabels = cellstr("LV Strain #" + string(1:5));
plot_BestLHS_Strains(lvPROBinfo,ABXinfo,LVStrain_select_outcomes_4wk,...
    LVStrain_all_warnall_4wk,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])

%% 9a. LV Strains every day
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ABXinfo.addABX = true; % true: administer ABX, false: no ABX
ABXinfo.dose = -2.64; % decay rate
ABXinfo.duration = 5; % days

% probiotic dosing regimen
wk1 = [1 1 1 1 1 1 1]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([1 1 1 1 1 1 1],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;
PROBinfo.numberDoses = 1; % days
PROBinfo.doseFrequency = NaN; % every XX days


PROBinfo.customEvalPoints = [-10:0.25:0.5, 1:0.25:6.5, 7:7:84, 77+30, 77+90, 77+180];
LVStrain_select_outcomes_7wk = NaN(size(LV_strains,1),...
    length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
LVStrain_all_warnall_7wk = NaN(size(LV_strains,1),length(POPinfo.Indexes));
for combo_id = 1:size(LV_strains,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    sel_strain = LV_strains(combo_id,:);

    PROBinfo.strainCharacteristics = sel_strain;
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    LVStrain_select_outcomes_7wk(combo_id,:,:,:) = select_outcomes;
    LVStrain_all_warnall_7wk(combo_id,:) = warnall == "";
end

save(strcat(result_fdr,'/F5/F5_LV_Best_Strains.mat'),...
    'LVStrain_all_warnall_7wk','LVStrain_select_outcomes_7wk','-append')


%% 9b. Plot Results everyday

% probiotic dosing regimen
wk1 = [1 1 1 1 1 1 1]; % wk 1: 1/day 4 days
wks = 10; % 2/wk for 10 wks
int_lbs = repmat([1 1 1 1 1 1 1],1,wks);
customDosing = find([wk1,int_lbs]); % get numerical values

lvPROBinfo = PROBinfo;
lvPROBinfo.customProbDosing = customDosing;

doseLabels = lvPROBinfo.customProbDosing+4;
definedStrain = LV_strains;

plotLabels = cellstr("LV Strain #" + string(1:5));
plot_BestLHS_Strains(lvPROBinfo,ABXinfo,LVStrain_select_outcomes_7wk,...
    LVStrain_all_warnall_7wk,definedStrain,doseLabels,plotLabels)
set(gcf,'units','inches','position',[0 0 16 10])







