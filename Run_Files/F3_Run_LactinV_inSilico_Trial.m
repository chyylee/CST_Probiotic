%% F3_Run_LactinV_inSilico_Trial.m
%
% GOAL: Replicate Lactin-V clinical trial and compare predicted BV
%   recurrences to clinically observed frequencies
%
% Study: Cohen et al. (2020) - Randomized Trial of Lactin-V to Prevent Recurrence of Bacterial Vaginosis
%
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. SET-UP RUN INFORMATION

% LOAD DEFAULT MODEL INPUTS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
clear; clc;
ws_name = 'input_files/SSConfig-Analysis-HMP-Virtual-Population.mat';
simulationType = "comboLV";
[POPinfo,PROBinfo,ABXinfo,~] = define_Common_Simulation_Inputs(ws_name,simulationType);
output_fdr_nm = 'result_workspaces/';


% DEFINE STRAINS TO EVAUATE
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

nullP = [0.5 0 0 0 0 0 0 -0.022];
modP = [0.5 -0.006  -0.04 -0.04 -0.012 -0.04 -0.04 -0.022];
badP = [0.5 -0.012  -0.04 -0.04 -0.006 -0.04 -0.04 -0.022];
tradP = [0.5 0 0 0 -0.012 0 0 -0.022];

strainList = {nullP,tradP,modP,badP};
run_nms = {'Null','Traditional','Moderate','Bad'};

%% 2. SET-UP CLINICAL REGIMENS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% CREATE ABX/PROBIOTIC REGIMENS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Requires input into structures for ABXinfo, PROBinfo, and POPinfo
%       - NOTE: POPinfo now has "SIMinfo" incorporated, rather than being
%               assigned it's own structure

ABXinfo.dose = -3.82;
LactinV = {ABXinfo,PROBinfo,POPinfo};

% ~~~~~~~~~~~~~~~~~ ABX + PLACEBO ~~~~~~~~~~~~~~~~~
PROBinfo.dose = 0; % relative unit

Placebo = {ABXinfo,PROBinfo,POPinfo};

% ~~~~~~~~~~~~~~~~~ COMPILE DESIGNS INTO CELL ARRAY ~~~~~~~~~~~~~~~~~
TestDesigns = {LactinV, Placebo};

%% 3. RUN TEST DESIGNS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ITERATE THROUGH EACH STRAIN AND REGIMEN TYPE
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

t = tic;
testDesignNames = {'Lactin-V','Placebo'};
numPatients = size(POPinfo.virtualPatientParameters,1);
for k = 1:length(strainList) % strain iteration
    all_select_outcomes = NaN(length(TestDesigns),numPatients,8+length(PROBinfo.customEvalPoints),5);
    all_warnall = NaN(length(TestDesigns),numPatients);
    probiotic_params = strainList{k};
    
    run_name = run_nms{k};
    for i = 1:length(TestDesigns) % Regigmen iteration
        tmp = TestDesigns{i};
        ABXinfo = tmp{1};
        PROBinfo = tmp{2};
        POPinfo = tmp{3};

        PROBinfo.strainCharacteristics = probiotic_params;
        [select_outcomes, ~, warnall, ...
            ~] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo);
        all_select_outcomes(i,:,:,:) = select_outcomes;
        all_warnall(i,:) = warnall == "";
    end
    if isfolder(output_fdr_nm)
        save(strcat(output_fdr_nm,'F3_',run_name,'_LactinVPlacebo.mat'),'all_warnall','all_select_outcomes',...
            'testDesignNames','TestDesigns','probiotic_params')
    else
        mkdir(output_fdr_nm)
        save(strcat(output_fdr_nm,run_name,'-LactinVPlacebo.mat'),'all_warnall','all_select_outcomes',...
            'testDesignNames','TestDesigns','probiotic_params')
    end
end

tend = toc(t);
disp(strcat(num2str(tend/60), " minutes"))
%% 4. VISUALIZE RESULTS
output_fdr = 'result_workspaces/F3/LV_Trial_A382/';
ws_list = strcat('F3_',{'Traditional_LactinVPlacebo.mat',...
    'Null_LactinVPlacebo.mat',...
    'Moderate_LactinVPlacebo.mat',...
    'Bad_LactinVPlacebo.mat'});

figure
plot_LactinVSimulationvsObserved(output_fdr,ws_list)
%%


