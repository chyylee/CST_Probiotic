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
load(ws_name,'strainListLHS')

%% 2. SET-UP CLINICAL REGIMENS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% CREATE ABX/PROBIOTIC REGIMEN
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Requires input into structures for ABXinfo, PROBinfo, and POPinfo
%       - NOTE: POPinfo now has "SIMinfo" incorporated, rather than being
%               assigned it's own structure

[param_names,SSnms,sp_cols,time_names] = get_naming_terms();

wk = [1 1 1 1 1 1 1]; % wk 1: 1/day 4 days
customDosing = find(wk); % get numerical values
PROBinfo.custom = true;
PROBinfo.customProbDosing = customDosing;

ABXinfo.addABX = false; % true: administer ABX, false: no ABX

%% 3. RUN SYSTEMATIC ANALYSIS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ITERATE THROUGH EACH STRAIN (computation ~10hrs)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

all_select_outcomes = NaN(size(strainListLHS,1),length(POPinfo.Indexes),8+length(PROBinfo.customEvalPoints),5);
all_warnall = NaN(size(strainListLHS,1),length(POPinfo.Indexes));
c = 0;
for combo_id = 1:size(strainListLHS,1)
    tic
    disp(strcat("Combination #", num2str(combo_id)))
    PROBinfo.strainCharacteristics = strainListLHS(combo_id,:);
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
    save('F4_Systematic_Acute_Probiotic.mat','all_warnall','all_select_outcomes','POPinfo','combo_id',...
        'strainListLHS', 'PROBinfo','ABXinfo')
end



     