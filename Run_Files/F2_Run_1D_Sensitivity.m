%% F2_Run_1D_Sensitivity.m
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

%% 2. RUN 1D SENSITIVITY ANALYSIS


% DEFINE THE 1D PARAMETER INFO TO ALTER 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
base_params = [0.5 0 0 0 0 0 0 -0.022]; % null probiotic as base
n = 21; % number of iterationss
param_alterkg = linspace(0,1,n); % growth rate parameter range
param_alteris = linspace(-0.1,0.1,n); % interspecies interaction parameter range
param_alterii = -logspace(log10(0.004),log10(0.04),n); % interspecies interaction parameter range
output_ws_nm = 'F2_1D_Sensitivity.mat';


% RUN SIMULATION FOR ALL ITERATIONS OF PARAMETER CHANGES (NOTE: This code takes 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
numVirtualPatients = length(POPinfo.virtualPatientComposition);
all_select_outcomes = NaN(length(base_params),length(param_alterkg),numVirtualPatients,8,5);
all_warnall = NaN(length(base_params),length(param_alterkg),numVirtualPatients);
t = tic;
for param_id = 1:length(SIMinfo.alteredParamIndex)
    if param_id == 1
        param_alter = param_alterkg;
    elseif param_id == length(SIMinfo.alteredParamIndex)
        param_alter = param_alterii;
    else
        param_alter = param_alteris;
    end
    for pert_id = 1:length(param_alter)
        SIMinfo.paramAlteration = base_params;
        SIMinfo.paramAlteration(param_id) = param_alter(pert_id);
        [select_outcomes,nonBVflag,tall,yall,warnall,...
           dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo);

        all_select_outcomes(param_id,pert_id,:,:,:) = select_outcomes;
        all_warnall(param_id,pert_id,:) = warnall == "";


        disp(strcat("Parameter #",num2str(param_id), " iterations #",num2str(pert_id)))
    end
    save(strcat(output_fdr,output_ws_nm),'all_select_outcomes','all_warnall',...
        'POPinfo','PROBinfo','SIMinfo','base_params',...
        'param_alterkg','param_alteris','param_id','pert_id',...
        'sp_cols','param_names','param_alterii')
    pause(60*3) % pause for 3 minutes
end
tend = toc(t);
disp(stcat(num2str(tend/60),'min'))
 
%% 3. PLOT RESULT

% PULL DATA FOR PLOTTING (load workspace if needed)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
output_fdr = 'result_workspaces/';
load(strcat(output_fdr,'F2_1D_Sensitivity.mat'))
time_id = 8; % time point to evaluate
params_to_view = SIMinfo.alteredParamIndex;
figFlag = false; % true: plot parameters separately; false: subplot

[param_names,SSnms,sp_cols,time_names] = get_naming_terms();

plot_1D_SensitivityAnalysis(all_select_outcomes,all_warnall,time_id,...
    params_to_view,param_alterkg,param_alteris,param_alterii,param_names,SSnms,sp_cols,figFlag)

%%
function plot_1D_SensitivityAnalysis(all_select_outcomes,all_warnall,time_id,...
    params_to_view,param_alterkg,param_alteris,param_alterii,param_names,SSnms,sp_cols,figFlag)

%%
    param_alter = param_alteris;
    % Loop through each parameter, plot 1D
    for param_id = 1:size(all_select_outcomes,1)
        allCounts = NaN(length(param_alter),4);
        allFrequencies = NaN(length(param_alter),4);
        for pert_id = 1:length(param_alter)
            Xwarn = squeeze(all_warnall(param_id,pert_id,:)) == 1;
            X = squeeze(all_select_outcomes(param_id,pert_id,Xwarn,time_id,2:end));
            CST_post = get_CST_post(X);
            Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
            Frequency = Counts ./ sum(Counts,2);
            allCounts(pert_id,:) = Counts;
            allFrequencies(pert_id,:) = Frequency;
        end
        
        if figFlag
            figure
        else
            subplot(2,4,param_id)
        end
        Y = allCounts;
        Y = bsxfun(@rdivide, Y, sum(Y,2))*100;
        bar(Y,'stacked')
        xticks(1:length(param_alter))
        if param_id == 1
            xticklabels(param_alterkg) % first parameter is kgrow-Probiotic (different units)
        elseif param_id == 8
            xticklabels(param_alterii)
        else
            xticklabels(param_alter)
        end 
        colororder(sp_cols)
        xlabel('Parameter Value')
        ylabel('Frequency')
        title(param_names(params_to_view(param_id)))
        xtickangle(90)
        set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
        ylim([0 100])
    end
    legend(SSnms)
end


