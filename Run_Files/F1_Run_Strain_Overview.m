%% F1_Strain_Overview.m
% 
% GOAL: Demonstrate simulation of probiotic dosing for a 'null' probiotic
% strain and a 'traditional' probiotic strain
%   Figure 1: Time series plot of composition before, during and after
%       probiotic
%   Figure 2: Bar chart of the response type evaluated at different time
%       points
%
% REQUIRES: 
%  - Model inputs as defined in define_Common_Simulation_Inputs
%       - Defines virtual patient population, simulation inputs
%  - simulate_CST_probiotic_object.m
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Author: Christina Lee
% v1: April 5th, 2023
% v2: April 26th, 2023 (updated to include P->P interactions as a probiotic 
%       characteristic)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1.LOAD MODEL DEFAULTS
clear; clc;
ws_name = 'input_files/SSConfig-Analysis-HMP-Virtual-Population.mat';
simulationType = "base";
[POPinfo,PROBinfo,~,SIMinfo] = define_Common_Simulation_Inputs(ws_name,simulationType);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% POPinfo: struct that defines virtual patient population
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   - virtualPatientParameters: # patients x 20 parameter matrix generated
%           by LHS
%   - virtualPatientComposition: 1 x # patients cell of analytically
%            predicted steady-state absolute abundances of nAB, Li, oLB
%            (predicted without considering probiotic in community)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PROBinfo: struct that defines probiotic parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   - dose: absolute abundance of probiotic when dosed into the base
%           community (selected from relative abundance data observed in
%            Dausset et al. 2018)
%   - numberDoses: total number of doses given
%   - doseFrequency: time separation (days) between each dose
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SIMinfo: struct that defines simulation inputs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   - startProbiotic: time point after simulation begins to start probiotic
%           therapy
%   - stopSimulation: time point to end simulation relative to the last
%           day of probiotic treatment
%   - initialComposition: scaling factor to "displace" from the predicted
%           composition (POPinfo.virtualPatientComposition). For example,
%           -0.01 is a 1% decrease in abundance and 0.05 is a 0.5% increase
%           from the predicted abundance
%   - alteredParamIndex: indexes of probiotic-related parameters relative
%           to the full parameter list 
%   - paramAlteration: new value for probiotic parameters (a 1 x length
%           alteredParamIndex array of values)
%   - paramAlterationType: 
%           'None' (use LHS probiotic parameter, ignore paramAlteration values),
%           'fold addition' (scale LHS probiotic parameter (LP) + paramAlteration * LP, 
%           'absolute' (set to exact value in paramAlteration)
%   - ODEoptions: "options" input for ode solvers
%   - ODEsolver: ODE solver function handle (e.g., @ode45, @ode15s)

%% 2. RUN NULL PROBIOTIC MODEL

SIMinfo.paramAlteration = [0.5 0 0 0 0 0 0 -0.0220]; % null probiotic (average values) - 2 ODE errors, 99 seconds

tic
[select_outcomes,nonBVflag,tall,yall,warnall,...
    dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo);
toc

% save results
output_fdr = 'result_workspaces/';
if ~isfolder(output_fdr)
    mkdir(output_fdr)
end
fnm = strcat(output_fdr,'F1_Null_Strain.mat');
save(fnm,"select_outcomes","nonBVflag","tall","yall","warnall","dose_lbls",...
    "POPinfo","PROBinfo","SIMinfo")

%% 3. RUN TRADITIONAL PROBIOTIC MODEL

SIMinfo.paramAlteration = [0.5 0 0 0 -0.01 0 0 -0.0220]; % null + strong inhibition of nAB 324.576967 seconds.

tic
[select_outcomes,nonBVflag,tall,yall,warnall,...
    dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo);
toc

fnm = strcat(output_fdr,'F1_Traditional_Strain.mat');
save(fnm,"select_outcomes","nonBVflag","tall","yall","warnall","dose_lbls",...
    "POPinfo","PROBinfo","SIMinfo")

%% 4. RUN NON-TRADITIONAL PROBIOTIC MODEL

SIMinfo.paramAlteration = [0.5 0.01 0 0 0 0 0 -0.0220]; % null + strong inhibition of P 324.576967 seconds.

tic
[select_outcomes,nonBVflag,tall,yall,warnall,...
    dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo);
toc

fnm = strcat(output_fdr,'F1_NonTraditional_Strain.mat');
save(fnm,"select_outcomes","nonBVflag","tall","yall","warnall","dose_lbls",...
    "POPinfo","PROBinfo","SIMinfo")

%% 4. PLOT RESULTS
output_fdr = 'result_workspaces/F1/';
ws_names = {'F1_Null_Strain.mat', 'F1_Traditional_Strain.mat','F1_NonTraditional_Strain.mat'};
[~,SSnms,sp_cols,time_names] = get_naming_terms();


null = load(strcat(output_fdr,ws_names{1}));
trad = load(strcat(output_fdr,ws_names{2}));
ntrad = load(strcat(output_fdr,ws_names{3}));

f1a(1) = figure;
time_id = 8; % 12 mo
for resp_id = 1:4
    subplot(2,2,resp_id)
    plot_ProbioticTrajectoryByResponse(resp_id,null.select_outcomes,...
        null.warnall,null.tall,null.yall,time_id,null.dose_lbls,SSnms,sp_cols);
end

f1(1) = figure;
selected_timepoints = [2 3 5:8]; % time points evaluated (0mo - 12mo post)
subplot(2,4,[1 5])
plot_ResponseFrequenciesByTimepoint(selected_timepoints,null.select_outcomes,...
    null.warnall,SSnms,time_names,sp_cols);
title(strrep(extractBetween(ws_names{1},'1_','.mat'),"_", " "))
legend(SSnms,'location','southoutside','Orientation','horizontal')
subplot(2,4,[2 6])
plot_ResponseFrequenciesByTimepoint(selected_timepoints,trad.select_outcomes,...
    trad.warnall,SSnms,time_names,sp_cols);
title(strrep(extractBetween(ws_names{2},'1_','.mat'),"_", " "))
legend(SSnms,'location','southoutside','Orientation','horizontal')
idx = [3 4 7 8];
for resp_id = 1:4
    subplot(2,4,idx(resp_id))
    plot_ProbioticTrajectoryByResponse(resp_id,trad.select_outcomes,...
        trad.warnall,trad.tall,trad.yall,time_id,trad.dose_lbls,SSnms,sp_cols);
end
set(gcf,'units','inches','position', [0 0 18 5])

f1a(2) = figure;
plot_ResponseFrequenciesByTimepoint(selected_timepoints,...
    ntrad.select_outcomes,ntrad.warnall,SSnms,time_names,sp_cols);
title(strrep(extractBetween(ws_names{3},'1_','.mat'),"_", " "))

