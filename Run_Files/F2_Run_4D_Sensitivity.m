%% F2_Run_4D_Sensitivity.m
%
% GOAL: Generate 1D sensitivity analysis on probiotic parameters starting with
% the "null" probiotic. The "null" probiotic is defined as having a
% moderate growth rate (0.5) and all interspecies interaction terms as
% zero.
%
% REQUIRES: 
%  - MATLAB workspace with:
%       - LHSmat, classification of LHSmat steady-state types
%       - Output will start with 'SSConfig-Analysis..."
%  - simulate_CST_probiotic_object.m
%
% OUTPUT:
%   - Workspace with 1D analysis results saved
%   - Figure plotting the response frequency for each parameter alteration
%       at a given evaluation time point
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1.LOAD MODEL DEFAULTS
clear; clc;
output_fdr = 'result_workspaces/';
ws_name = 'input_files/SSConfig-Analysis-HMP-Virtual-Population.mat';
simulationType = "base";
[POPinfo,PROBinfo,~,SIMinfo] = define_Common_Simulation_Inputs(ws_name,simulationType);
[param_names,~,sp_cols,~] = get_naming_terms();
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
%% 2. RUN 4D SENSITIVITY ANALYSIS


% DEFINE THE 4 PARAMETERS TO ALTER
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
base_params = [0.5 0 0 0 0 0 0 -0.022]; % base parameter set
pmags1 = [-0.01 0 0.01]; % magnitudues for the parameter alteration
pidx = [8 17 18 19]; % indexes of parameter combinations
output_ws_nm = 'F2_4D_Sensitivity.mat';

vectorCell = {[pmags1],[pmags1],[pmags1],[pmags1]};
[newValueMat] = generate_input_combos(vectorCell);
fullparamChange = zeros(size(newValueMat,1),length(base_params));
fullparamChange(:,1) = base_params(1);
fullparamChange(:,end) = base_params(end);
fullparamChange(:,[2,5,6,7]) = newValueMat;

% RUN ALL COMBINATIONS OF ALTERATION
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
numVirtualPatients = length(POPinfo.virtualPatientComposition);
all_select_outcomes = NaN(size(newValueMat,1),numVirtualPatients,8,5);
all_warnall = NaN(size(newValueMat,1),numVirtualPatients);
c = 1;
for pert_id = 1:size(newValueMat,1)
    SIMinfo.paramAlteration = fullparamChange(pert_id,:);
    [select_outcomes,~,~,~,warnall,~] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo);
    all_select_outcomes(pert_id,:,:,:) = select_outcomes;
    all_warnall(pert_id,:) = warnall == "";


    % Pause simulations for 1 minute, every 10 iterations
    c = c + 1;
    if c == 10
        save(strcat(output_fdr,output_ws_nm),'all_select_outcomes','all_warnall',...
            'POPinfo','PROBinfo','SIMinfo','base_params',...
            "all_warnall", 'pert_id','newValueMat','pidx')
        c = 1; pause(60);
    end
end

%% 3. PLOT RESULTS
output_fdr = 'result_workspaces/';

load(strcat(output_fdr,'F2_4D_Sensitivity_Apr28.mat'))
time_id = 8;

[allCounts,pChiDfAllResponses, pChiDfBvResponses] = ...
    plot_4D_SensitivityAnalysis(all_select_outcomes,all_warnall,...
    newValueMat,pidx,time_id,param_names,sp_cols)
ylim([0,100])



