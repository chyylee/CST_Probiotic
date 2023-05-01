%% F2_Run_Best_Designs.m
%
% GOAL: From 1D results, pick what would have the best efficacy with the smallest
% perturbation
%
%   INCREASING EFFICACY
%           nAB -> P, +0.012
%           P -> nAB, -0.012
%           nAB -> P, +0.012 & P -> nAB, -0.012
%           nAB -> P, +0.012 & P -> nAB, -0.012 & kgrow-P = 0.8
%      Rev: nAB -> P, -0.012 & P -> nAB, -0.012
%
%
%   BOOSTING oLB over Li/Probiotic
%       P -> oLB, +0.012
%       P -> oLB, +0.012 & P -> nAB, -0.012
%       P -> oLB, +0.012 & P -> nAB, -0.012 & P -> Li, -0.012
%
% REQUIRES: 
%  - MATLAB workspace with:
%       - LHSmat, classification of LHSmat steady-state types
%       - Output will start with 'SSConfig-Analysis..."
%  - simulate_CST_probiotic_object.m
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

%% 2. RUN "DESIGNED" PROBIOTICS

pidx_list = {1,2, 5, [2 5], [1 2 5], [2 5], 7, [5 7],[5 6 7]}; % parameter indexes
mag = 0.01; % magnitude
kgr = 0.8; % final growth rate
pert_list = {NaN,1,-1,[1 -1],[kgr/mag 1 -1],[-1 -1],1,[-1 1],[-1 -1 1]}; % directionality
base_params = [0.5 0 0 0 0 0 0 -0.022]; % Null Probiotic
output_ws_nm = 'F2_Select_Probiotic_Designs.mat';

numVirtualPatients = length(POPinfo.virtualPatientComposition);
all_select_outcomes = NaN(length(pidx_list),numVirtualPatients,8,5);
all_warnall = NaN(length(pidx_list),numVirtualPatients);
for pert_id = 1:length(pidx_list)
    SIMinfo.paramAlteration = base_params; 
    if ~isnan(pert_list{pert_id})
        SIMinfo.paramAlteration(pidx_list{pert_id}) = mag*pert_list{pert_id};
    end

    [select_outcomes,nonBVflag,tall,yall,warnall,...
           dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo);

    all_select_outcomes(pert_id,:,:,:) = select_outcomes;
    all_warnall(pert_id,:) = warnall == "";

    save(output_ws_nm,'all_select_outcomes','all_warnall',...
        'POPinfo','PROBinfo','SIMinfo','base_params',...
        'pidx_list','pert_list','mag','kgr','pert_id')
end

%% 3. PLOT RESULTS
output_fdr = 'result_workspaces/';

load(strcat(output_fdr,'F2_Select_Probiotic_Designs_Apr28.mat'))
time_id = 8; % timepoint (8 = 12mo)
gr1 = [0 2 1 3 4 5]+1; % indexes for first plot
gr2 = [0 6:8]+1; % indexes for second plot
pert_names = {'null','1p: Best','1p: Traditional','2p: Best','3p: Best','2p: Counter','1p: oLB','2p: oLB','3p: oLB'};
    

[allCounts,pChiDfAllResponses, pChiDfBvResponses] = plot_SelectProbioticDesigns(pidx_list,...
    pert_list,pert_names,gr1,gr2,all_select_outcomes,all_warnall,...
    time_id);

