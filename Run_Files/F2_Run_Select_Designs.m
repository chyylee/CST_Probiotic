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

%%
function [allCounts,pChiDfAllResponses, pChiDfBvResponses] = plot_SelectProbioticDesigns(pidx_list,...
    pert_list,pert_names,gr1,gr2,all_select_outcomes,all_warnall,time_id)

    [param_names,~,sp_cols,~] = get_naming_terms();
    
    % PULL DATA FOR PLOTTING 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    allCounts = NaN(length(pidx_list),4);
    allFrequencies = NaN(length(pidx_list),4);
    param_vals = zeros(7,length(pidx_list));
    for pert_id = 1:length(pidx_list)
        abundance = squeeze(all_select_outcomes(pert_id,all_warnall(pert_id,:)==1,time_id,2:5));
        CST_post = get_CST_post(abundance);
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts,2);
        allCounts(pert_id,:) = Counts;
        allFrequencies(pert_id,:) = Frequency;
        
        param_vals(pidx_list{pert_id},pert_id) = pert_list{pert_id};
    end
    
   
    % STATISTICS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pChiDfAllResponses = NaN(length(pidx_list),length(pidx_list),3);
    pChiDfBvResponses = NaN(length(pidx_list),length(pidx_list),3);
    for pert_id1 = 1:length(pidx_list)
        for pert_id2 = 1:length(pidx_list)
            observed = [allCounts(pert_id1,:)', allCounts(pert_id2,:)'];
            [p, chi2stat,df] = chigof(observed);
            pChiDfAllResponses(pert_id1,pert_id2,:) = [p, chi2stat,df]; % Across all 4 reseponse types
        
            observed_BV = [observed(1,:);sum(observed(2:end,:))];
            [p_BV, chi2stat_BV,df_BV] = chigof(observed_BV);
            pChiDfBvResponses(pert_id1,pert_id2,:) = [p_BV, chi2stat_BV,df_BV]; % BV+/-
        end
    end
    
    
    % PLOT THE RESULTS 
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    % ~~~~~~~~~~~~ Plots Parameter Alterations  ~~~~~~~~~~~~ 
    subplot(4,2,7)
    param_vals(param_vals > 1) = 1;
    ax = imagesc(param_vals([1,2,5],gr1));
    colormap(gca,[0 0.75 0.75; ...
        1 1 1; ...
        0.75 0 0.75])
    set(gca,'YTICK',1:3,'YTICKLABELS',param_names([4,8,17]),...
        'XTICK',1:length(gr1),'XTICKLABELS',pert_names(gr1))
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    xtickangle(90)
    
    subplot(4,2,8)
    param_vals(param_vals > 1) = 1;
    ax = imagesc(param_vals([5,6,7],gr2));
    colormap(gca,[0 0.75 0.75; ...
        1 1 1; ...
        0.75 0 0.75])
    set(gca,'YTICK',1:3,'YTICKLABELS',param_names([17,18,19]),...
        'XTICK',1:length(gr2),'XTICKLABELS',pert_names(gr2))
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    xtickangle(90)
    
    % ~~~~~~~~~~~~ Plots the BV Design Comparisons  ~~~~~~~~~~~~ 
    subplot(4,2,[1,3,5])
    X = allFrequencies(gr1,:)*100;
    bar(X,'stacked')
    colororder(sp_cols)
    xlim([0.5 length(gr1)+0.5])
    xticklabels([])
    ylabel('Percent Samples')
    hold on
    Y = ones(length(gr1),1)*100;
    P = squeeze(pChiDfBvResponses(gr1,gr1,1));
    superbar(Y, 'P',round(P,4),'BarFaceColor', 'none', 'BarEdgeColor', 'none')
    ylim([0 200])
    
    for i = 1:size(X,1) % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f %%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center',...
            'fontsize',10)
    end
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    title('Reducing BV+')
    
    % ~~~~~~~~~~~~ Plots the LB Design Comparisons  ~~~~~~~~~~~~ 
    subplot(4,2,[2,4,6])
    X = allFrequencies(gr2,:)*100;
    bar(X,'stacked')
    colororder(sp_cols)
    xlim([0.5 length(gr2)+0.5])
    xticklabels([])
    ylabel('Percent Samples')
    hold on
    Y = ones(length(gr2),1)*100;
    P = squeeze(pChiDfBvResponses(gr2,gr2,1));
    superbar(Y, 'P',round(P,4),'BarFaceColor', 'none', 'BarEdgeColor', 'none')
    ylim([0 200])
    
    for i = 1:size(X,1) % Label for percent nAB response post
        text(i,X(i,1),num2str(X(i,1),'%.1f %%'),'VerticalAlignment','top', 'HorizontalAlignment', 'center',...
            'fontsize',10)
    end
    set(gca,'XColor','k','YColor','k','LineWidth',1,'fontsize',10)
    title('Boosting oLB')
end