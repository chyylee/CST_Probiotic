%% F4_Run_Systematic_LactinV_Probiotic.m
%
% GOAL: Methodically test te effect of a variety of different probiotic 
% "strains" across an in silico population (matched for recurrent BV)
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. SET-UP RUN INFORMATION

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% LOAD DEFAULT MODEL INPUTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc;
output_fdr = 'result_workspaces/';
ws_name = 'input_workspaces/SSConfig-Analysis-Virtual-Population.mat';
simulationType = "comboLV";
[POPinfo,PROBinfo,ABXinfo,~] = define_Common_Simulation_Inputs(ws_name,simulationType);
output_fdr_nm = 'result_workspaces/';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% DEFINE STRAINS TO EVAUATE
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% newValueMat (previously defined by LHS)
%   NOTE: newValueMat is a # probiotic strains x 7 matrix  generated by LHS
load(ws_name,'newValueMat')


%% 3. RUN SYSTEMATIC ANALYSIS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ITERATE THROUGH EACH STRAIN (computation ~10hrs)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_select_outcomes = NaN(size(newValueMat,1),length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
all_warnall = NaN(size(newValueMat,1),length(POPinfo.Indexes));
c = 0;
for combo_id = 1:size(newValueMat,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    PROBinfo.strainCharacteristics = newValueMat(combo_id,:);
    [select_outcomes, meanstdDuring, warnall, ...
        nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);

    all_select_outcomes(combo_id,:,:,:) = select_outcomes;
    all_warnall(combo_id,:) = warnall == "";

    % After ten runs, pause for 5 minutes
    disp(strcat("Combination #", num2str(combo_id)))
    if c == 10
        pause(60*5); c = 0;
    else
        c = c + 1;
    end
    toc

    % Saves results with each iteration (code takes ~10 hours to run,
    % protects against losing data)
    save('F4_Systematic_LactinV_Probiotic.mat','all_warnall','all_select_outcomes','POPinfo','combo_id',...
        'newValueMat', 'PROBinfo','ABXinfo')
end

%% 4. Load workspace generated in Fig4_Systematic_LactinV_Probiotic.m
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PROCESS RESULTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_fdr_nm = 'result_workspaces/F4/';
load(strcat(output_fdr_nm,'F4_Systematic_LactinV_Probiotic.mat'))
    % all_select_outcomes: Numbers Strains x Simulated Patients x Time
    %       Points x [time point, nAB, Li, oLB, Probiotic] Abundance
    % all_warnall: Indicates runs without errors (1 = no error, 0 = error)
    % strainListLHS: Probiotic strain parameter values (Numbers Strains x
    %       Strain Parameters)
    % POPinfo: Structure with information on simulated population
    % PROBinfo: Structure with information about probiotic
 

%% 2. Pull Data of Strains at Lactin-V Evaluation Points
% NOTE: Lactin-V is evaluated at 12 weeks and 24 weeks after initial
% probiotic dosing. This corresponds to index 12 and 14 of
% all_select_outcomes

[param_names,SSnms,sp_cols,time_names] = get_naming_terms();
time_id = [12 14]; % 12 or 24 weeks
numPatients = size(all_select_outcomes,2);
numStrains = size(strainListLHS,1);
numODEwarnings = sum(all_warnall==0,2);

errorThresh = 1; % must have a least one successful run (no ODE solver errors)

% collect information about response counts and frequencies
strainCounts12 = NaN(numStrains,4);
strainFrequency12 = NaN(numStrains,4);
strainCounts24 = NaN(numStrains,4);
strainFrequency24 = NaN(numStrains,4);

% loop through each "strain" (strainListLHS combination)
for strain_id = 1:numStrains
    warn_id = all_warnall(strain_id,:) == 1;
    if sum(warn_id) > errorThresh
        CST_post = get_CST_post(squeeze(all_select_outcomes(strain_id,warn_id,time_id(1),2:end)));
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts);
        strainCounts12(strain_id,:) = Counts;
        strainFrequency12(strain_id,:) = Frequency;

        CST_post = get_CST_post(squeeze(all_select_outcomes(strain_id,warn_id,time_id(2),2:end)));
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts);
        strainCounts24(strain_id,:) = Counts;
        strainFrequency24(strain_id,:) = Frequency;
    end
end
disp('Maximum Response Frequency 12wks: nAB, Li, oLB, Probiotic:')
disp(max(strainFrequency12))

disp('Maximum Response Frequency 12wks: nAB, Li, oLB, Probiotic:')
disp(max(strainFrequency24))

%% 3. DETERMINE WHICH PROBIOTICS MATCH LACTIN-V RESULTS
ClinicalRecur = [46 59]; % Number of patients that were BV+ by week 12, week 24
ClinicalCure = [87 63]; % Number of patients that were not BV+ by week 12, week 24
ClinicalN = ClinicalCure + ClinicalRecur; % Total number of patients with known outcome
ClinicalResults = (ClinicalRecur./ClinicalN)*100; % Frequency of Recurrence (BV+/nAB dominance)

% Loop through results, collect distance from clinical result
tot_dist = NaN(size(strainFrequency12,1),1);
for i = 1:size(strainFrequency12,1)
    tmp12 = strainFrequency12(i,1)*100; % Frequency of nAB dominated states at 12 weeks
    tmp24 = strainFrequency24(i,1)*100; % Frequency of nAB dominated states at 24 weeks

    d12 = abs(tmp12-ClinicalResults(1)); % absolute distance at 12 week
    d24 = abs(tmp24-ClinicalResults(2)); % absolute distance at 24 week
    
    tot_dist(i) = d12 + d24; % total distance
end

errorThresh = 500; % threshold for the number of successful runs (244 gives 1000 "good" runs)

er_id = numODEwarnings <= errorThresh;
histogram(tot_dist(er_id))

% quick viz of results:
title('Absolute Distance from Clinical Observation')
ylabel('Frequency')

% select acceptable deviation from clinical results
dis_id = tot_dist < 5; % total percent different
LV_strains = strainListLHS(dis_id & er_id,:);

save(strcat(output_fdr_nm,'F4_Systematic_LactinV_Probiotic.mat'),'LV_strains','tot_dist','-append')

%% 4. PLOT CANDIDATE LV STRAIN PARAMETERS
subplot(1,4,1)
boxchart(LV_strains(:,1));
hold on
x = repmat(1,size(LV_strains,1),1);
s1 = swarmchart(x,LV_strains(:,1),[],'red','XJitterWidth',0.5);
xticklabels(param_names(PROBinfo.strainCharIndexes(1)))
ylabel('Parameter Value')
title('Growth Term')

subplot(1,4,[2,4])
boxchart(LV_strains(:,2:end))
hold on
x = repmat(1:size(LV_strains,2)-1,size(LV_strains,1),1);
s2 = swarmchart(x,LV_strains(:,2:end),[],'red','XJitterWidth',0.5);
xticklabels(param_names(PROBinfo.strainCharIndexes(2:end)))
ylabel('Parameter Value')
hold on
yline(0,'k')
title('Interaction Terms')

set(gcf,'units','pixels','position',[1 1 900 300])


%% Visualize Top Performing Strains
X = strainFrequency24; % or strainFrequency12

disp('Top Performing Strain - Percentatge of Subjects with nAB, Li, oLB, and P-dominance')
disp(prctile(X(er_id,:)*100,100))

sortThresh = prctile(X(er_id,:),90); % analyze 90th percentile strains

[~,Response_Label] = max(X,[],2);
Response_Label(isnan(X(:,1))) = -1;
bestStrains = cell(5,1);
bestStrains(5) = {strainListLHS(tot_dist == min(tot_dist(er_id)),:)};
for resp_id = 1:4
    ridx =  Response_Label == resp_id;
    eidx = numODEwarnings <= errorThresh;
    aidx = X(:,resp_id) > sortThresh(resp_id);

    selected_strains = strainListLHS(ridx & eidx & aidx,:);

    Xtmp = max(X(ridx & eidx & aidx,resp_id));

    bestStrains(resp_id) = {strainListLHS(ridx & eidx & X(:,resp_id) == Xtmp,:)};

    figure;
    subplot(1,4,1)
    boxchart(selected_strains(:,1));
    hold on
    x = repmat(1,size(selected_strains,1),1);
    s1 = swarmchart(x,selected_strains(:,1),[],'red','XJitterWidth',0.5);
    xticklabels(param_names(PROBinfo.strainCharIndexes(1)))
    ylabel('Parameter Value')
    title('Growth Term')
    
    subplot(1,4,[2,4])
    boxchart(selected_strains(:,2:end))
    hold on
    x = repmat(1:size(selected_strains,2)-1,size(selected_strains,1),1);
    s2 = swarmchart(x,selected_strains(:,2:end),[],'red','XJitterWidth',0.5);
    xticklabels(param_names(PROBinfo.strainCharIndexes(2:end)))
    ylabel('Parameter Value')
    hold on
    yline(0,'k')
    title('Interaction Terms')
    
    set(gcf,'units','pixels','position',[1 1 900 300])
end

save(strcat(output_fdr_nm,'F4_Systematic_LactinV_Probiotic.mat'),'bestStrains','-append')

%% RUN PLS-DA ANALYSIS
errorThresh = 500;
X = strainFrequency24; % or strainFrequency12

[~,lb] = max(X,[],2); % classify by nAB/Li/oLB/P is dominant response profile
lb(isnan(X(:,1))) = -1; % remove entries with NaN

prism_dat = [];
Xblock = [];
yblock = [];
c = 1;
for i = [1 2 3 4]
    tmp = strainListLHS(lb == i & numODEwarnings <= errorThresh,:);
    Xblock = [Xblock;tmp];
    yblock = [yblock;repmat(c,size(tmp,1),1)];
    c = c + 1;
end

% ~~~~~~~~ Uncomment to Compare BV+ / BV-
% yblock(yblock == 3) = 2;
% yblock(yblock == 4) = 2;
% classes = {'nAB','LB'};
% filettl = '___________-24wk-LV-BVefficacy';
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

classes = {'nAB','Li','oLB','P'};
ynames = param_names(PROBinfo.strainCharIndexes);
Yblock = dummyvar(yblock);
filettl = '___________-24wk-LV-ResponseTypes';

ttl = '';
yname = '';
ep = 1;
[plsda_model,cls_error,cls_error_CV,ellip_points] = fPLSDA(Xblock,...
    Yblock,ynames,filettl,classes,yname,ttl,1,ep);

%% Plot PLSDA Results
output_fdr = 'result_workspaces/F4/LV_PLSDA/';
ws_list = {'30-Apr-2023-24wk-LV-BVefficacy_PLSDA_Model.mat',...
    '30-Apr-2023-24wk-LV-ResponseTypes_PLSDA_Model.mat', ...
    '30-Apr-2023-12wk-LV-BVefficacy_PLSDA_Model.mat',...
    '30-Apr-2023-12wk-LV-ResponseTypes_PLSDA_Model.mat'};

[~,~,sp_cols,~] = get_naming_terms();
for i = 1:length(ws_list)
    load(strcat(output_fdr,ws_list{i}))

    ttl = strrep(extractBetween(ws_list{i},'2023-','_PLSDA'),'-'," ");
    filettl = '';
    
    if length(classes) == 2
        colors = [brighten(sp_cols(1,:),-0.5); 0 0.4470 0.7410];
    else
        colors = sp_cols;
    end

    plot_PLS(LV1_scores, LV2_scores, LV1_loads,LV2_loads,...
            yblock, scores_xlabel, scores_ylabel, loads_xlabel, loads_ylabel,...
            classes,  vip_xnames, "", ttl, filettl,...
            errortext,0, scores_num_LVs, ellip_points,colors);
end