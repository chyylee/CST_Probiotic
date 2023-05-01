%% F3_Run_ABX_Probiotic_Acute.m
%
% GOAL: Compare acute probiotic to antibiotics, or combination antibiotic
% plus acute probiotic
%
%
% REQUIRES: 
%  - MATLAB workspace with:
%       - LHSmat, classification of LHSmat steady-state types
%       - Output will start with 'SSConfig-Analysis..."
%  - simulate_Clinical_Regimens.m
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% 1. SET-UP RUN INFORMATION
clear; clc;
ws_name = 'input_files/SSConfig-Analysis-HMP-Virtual-Population.mat';
simulationType = "combo";
[POPinfo,PROBinfo,ABXinfo,~] = define_Common_Simulation_Inputs(ws_name,simulationType);
output_fdr = 'result_workspaces/';

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% POPinfo: struct that defines virtual patient population
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   - virtualPatientParameters: # patients x 20 parameter matrix generated
%           by LHS
%   - virtualPatientComposition: 1 x # patients cell of analytically
%            predicted steady-state absolute abundances of nAB, Li, oLB
%            (predicted without considering probiotic in community)
%   - Indexes: indexes of the virtualPatientPopulaion relative to the
%           initial LHS matrix defined in SSConfig-Analysis_XXXX.mat
%   - InitialRelativeAbundance: scaling factor to "displace" from the predicted
%           composition (POPinfo.virtualPatientComposition). For example,
%           -0.01 is a 1% decrease in abundance and 0.05 is a 0.5% increase
%   - odeOptions: "options" input for ode solvers
%   - odeSolver: ODE solver function handle (e.g., @ode45, @ode15s)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% PROBinfo: struct that defines probiotic parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   - dose: absolute abundance of probiotic when dosed into the base
%           community (selected from relative abundance data observed in
%            Dausset et al. 2018)
%   - numberDoses: total number of doses given
%   - doseFrequency: time separation (days) between each dose
%   - StrainCharIndexes: indexes of probiotic-related parameters relative
%           to the full parameter list 
%   - StrainCharacteristics: values for probiotic parameters
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ABXinfo: struct that defines simulation inputs
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   - AddABX: true if ABX is to be simulated, false if no ABX simulated
%   - Dose: decay rate of nAB (from Mayer et al. 2015)
%   - Duration: number of days for treatment


%% 2. SET-UP CLINICAL REGIMENS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% CREATE ABX/PROBIOTIC REGIMENS
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%   Requires input into structures for ABXinfo, PROBinfo, and POPinfo
%       - NOTE: POPinfo now has "SIMinfo" incorporated, rather than being
%               assigned it's own structure

% ~~~~~~~~~~~~~~~~~ ABX + ACUTE PROBIOTIC ~~~~~~~~~~~~~~~~~
ABXinfo.addABX = true; % true = add antibiotic, false = skip
ABXinfo.dose = -3.82; % decay rate
ABXinfo.duration = 5; % days

PROBinfo.dose = 7.5; % relative unit
PROBinfo.numberDoses = 7; % days
PROBinfo.doseFrequency = 1; % every XX days
PROBinfo.strainCharIndexes =  [4,8,12,16,17,18,19,20];
PROBinfo.strainCharacteristics = [0.5 0 0 0 0 0 0 -0.022];

ABXNullDesign = {ABXinfo,PROBinfo,POPinfo};

% ~~~~~~~~~~~~~~~~~ ABX + NO PROBIOTIC/PLACEBO ~~~~~~~~~~~~~~~~~
% all other simulation info is the same, only need to change dose to zero
PROBinfo.dose = 0; % relative unit 

ABXOnlyDesign = {ABXinfo,PROBinfo,POPinfo};

% ~~~~~~~~~~~~~~~~~ NO ABX + ACUTE PROBIOTIC ~~~~~~~~~~~~~~~~~
ABXinfo.addABX = false; % true = add antibiotic, false = skip
ABXinfo.dose = 0; % decay rate
ABXinfo.duration = 0; % days

PROBinfo.dose = 7.5; % relative unit
PROBinfo.numberDoses = 7; % days
PROBinfo.doseFrequency = 1; % every XX days
PROBinfo.strainCharIndexes =  [4,8,12,16,17,18,19, 20];
PROBinfo.strainCharacteristics = [0.5 0 0 0 0 0 0 -0.022];

NullOnlyDesign = {ABXinfo,PROBinfo,POPinfo};

% ~~~~~~~~~~~~~~~~~ COMPILE DESIGNS INTO CELL ARRAY ~~~~~~~~~~~~~~~~~
testDesigns = {NullOnlyDesign,ABXNullDesign,ABXOnlyDesign};

%% 3. RUN TEST DESIGNS

% ITERATE THROUGH EACH REGIMEN TYPE
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

output_ws_nm = 'F3_ABX_Probiotic_ShortTerm.mat';
testDesignNames = {'Only Probiotic','ABX + Probiotic','ABX + Placebo'};

t = tic;
numPatients = size(POPinfo.virtualPatientParameters,1);
all_select_outcomes = NaN(length(testDesigns),numPatients,8,5);
all_warnall = NaN(length(testDesigns),numPatients);
for i = 1:length(testDesigns)
    tmp = testDesigns{i};
    ABXinfo = tmp{1};
    PROBinfo = tmp{2};
    POPinfo = tmp{3};
    [select_outcomes, ~, warnall, ...
        ~] = simulate_Clinical_Regimens(ABXinfo,PROBinfo,POPinfo);

    all_select_outcomes(i,:,:,:) = select_outcomes;
    all_warnall(i,:) = warnall == "";
end
tend = toc(t);
disp(strcat(num2str(tend/60), " minutes"))

save(strcat(output_fdr,output_ws_nm), 'all_select_outcomes', 'all_warnall',...
    "testDesignNames",'testDesigns')


%% 4. VISUALIZE RESULTS
output_fdr = 'result_workspaces/';
load(strcat(output_fdr,'F3_ABX_Probiotic_ShortTerm.mat'))

timePoints = [2 5 6 7 8]; % timepoints to analyze
[Counts,all_PvaluesLB,all_PvaluesBV] = ...
    plot_ProbioticABXComparisons(all_select_outcomes,all_warnall,...
    testDesigns,testDesignNames,timePoints);

