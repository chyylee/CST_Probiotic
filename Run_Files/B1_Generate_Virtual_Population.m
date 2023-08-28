%% G1_Generate_Reference_Virtual_Patient.m
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Goal: Use the initial/base sensitivity analysis to generate a larger
% reference population that can be used to match population-level
% observations in equilibrium behaviors (EB).
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Description of code workflow:
%   (1) Load the base/initial workspace
%           - Input: 'SSConfig-Analysis-XXXX.mat'
%   (2) Use parameter sets of known EB to generate new parameter sets
%           - Input: LHSmat selected by EB type
%           - Output: Larger parameter matrix with ("Reference"
%           population")
%   (3) Match the clinically observed EB Frequencies
%           - Input: EB frequencies for HMP or Gajer Cohort (use result of
%               G1_Analyze_Clinical_Equilibrium_Behaviors.m)
%           - Output: New SSConfig-Analysis.mat that has the matched
%           population data [WARNING: may overwrite previous workpaces]
%
% REQUIRES: [Add folders to MATLAB path]
%       * analyze_global_SS.m
%       * calc_SS_stability.m
%           - requires "S" which is the models symbolic equations, "Jmat" which
%           is the Jacobian of the system of ODEs
%       * analyze_Global_CST_SS.m
%   * MATLAB TOOLBOXES: Parallel Computing Toolbox, Bioinformatics Toolbox,
%           Symbolic Math Toolbox
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee (chyylee@umich.edu)
% University of Michigan
% Apr 28, 2023
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%% 1. Generate Base Virtual Population (no Probiotic)
rng(1);
input_fdr = 'input_files';
flnm = strcat(input_fdr,'lhs_settings_input.xlsx'); % Excel file with necessary parameters
[params,~] = pull_LHS_parameters(flnm);

NR = 30000; % Number of "samples", our work uses 5000

paramMatrix = defineMatrix(params, NR,'parameter'); % Generates LHS parameter sets

analyze_Global_SS(paramMatrix,{'nAB','Li','oLB'}) % Call function that determines SS configurations

%% 2. MATCH THE POPULATION TO CLINICAL FREQUENCIES
fdr_loc = "input_files/";
load(strcat(fdr_loc,'SSConfig-Analysis-Base.mat'))

rng(2);
NumPats = 2000; % Total size of the simulated populatiion

% ~~~~~~~~~~ HMP EB Frequencies ~~~~~~~~~~ (uncomment to run HMP)
hmp101 = [0.60 0.31 0.09]; % actual HMP distribution

counts_per_EB = ceil(NumPats.*hmp101); % 1SS nAB, 1SS Li, 1SS oLB, 2SS Li or nAB, 2SS Li or oLB, 2SS oLB or nAB
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sel_cst_names = {'1SS: [NO] CST-IV';
    '2SS: [NO] CST-IV or [Li] CST-III';
    '2SS: [NO] CST-IV or [oLB] CST-I/II/V'};

sp_combos = [1 0 0;
    1 1 0;
    1 0 1];

% Define Relative Abunance Range for Dominant Species
ub = 0.999;
lb = 0.500;

% Remove set with negative growth rates:
neg_gr = sum(LHSmat(:,1:3) < 0,2) > 0;
pos_si = sum(LHSmat(:,[4 8 12]) > 0,2) > 0;

s = RandStream('mlfg6331_64'); 
c = 1;
final_population = NaN(NumPats,length(param_names));
for cst_id = 1:length(sel_cst_names)
    indx = find(all_nm_CST == sel_cst_names{cst_id} & ~neg_gr & ~pos_si);

    % Checking Relative Abundance of Sample
    tmp = StbleSS(indx);
    abundance = cell2mat([tmp{:}]');
    rel_abundance = abundance ./ sum(abundance,2);
    ssu = find(sp_combos(cst_id,:));
    spComp = NaN(size(tmp,1),length(ssu));
    ssn = repmat(ssu,size(tmp,1),1);
    for i = 1:length(ssu)
        rtmp = rel_abundance(ssn==ssu(i),ssu(i));
        bid = rtmp > lb & rtmp < ub;
        spComp(:,i) = bid;
    end
    poss_indx = sum(spComp,2) == 1;
    % End Checking Relative Abundance of Sample

    mix_idx = indx(poss_indx); % pull idxes of all sets that meet ub/lb
    leftover_indx = indx(~poss_indx); % idxes of sets that do not meet ub/lb
    
    num_samples = counts_per_EB(cst_id); % number of samples needed given input EB distribution
    leftover_num_samples = num_samples - length(mix_idx); % Determine number of samples needed to be selected that do not meet ub/lb
    
    if leftover_num_samples >= 0
        tmp_id = randsample(s,leftover_indx,leftover_num_samples); 
    else
        tmp_id = randsample(s,mix_idx,num_samples); % if number of ub/lb sets is greater than the number of sets needed, randomly select
        mix_idx = [];
    end

    final_population(c:c+num_samples-1,:) = LHSmat([mix_idx;tmp_id],:); % final virtual population
    c = c + num_samples;
end

% Check Distribution Matches Desired Output (Rename
% "SSConfig-Analysis-Virtual-Population.mat"
analyze_Global_SS(final_population,{'nAB','Li','oLB'}) % Call function that determines SS configurations

%% 6. Define Probiotic Strains
rng(3);
input_fdr = 'input_files/';
flnm = strcat(input_fdr,'lhs_settings_input_probiotic.xlsx'); % Excel file with necessary parameters
[params,ICs] = pull_LHS_parameters(flnm);

NR = 500; % Number of "samples", our work uses 1000

pidx = [4,8,12,16,17,18,19,20]; % only pull the parameters associated with probiotic
paramMatrix = defineMatrix(params, NR,'parameter');
newValueMat = paramMatrix(:,pidx);
newValueMat(:,end) = repmat(-0.022,NR,1);

strainListLHS = newValueMat;

save(strcat(input_fdr,'SSConfig-Analysis-LV-Virtual-Population.mat'),'strainListLHS','-append')
