
function [POPinfo,PROBinfo,ABXinfo,SIMinfo] = define_Common_Simulation_Inputs(ws_name,simulationType)
   
    % SIMULATION INPUTS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    displace_SS = [-0.01 0.005 0.005 0]; % Change from predicted composition at SS (Initial Conditions)
    odeSolver = @ode15s; % ODE Solver
    odeOptions = odeset('AbsTol',1e-25,'RelTol',1e-10,'Jacobian',@gLV_jacobian); % Special ODE Solver Instructions
    null_probiotic = [0.5 0 0 0 0 0 0 -0.022]; % parameters values for null probiotic
    probiotic_param_indexes = [4,8,12,16,17,18,19,20]; % indexes in full parameter list for null probiotic
    
   
    % DEFAULT VIRTUAL PATIENT POPULATION
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    load(ws_name,'LHSmat','StbleSS')

    sel_nets_base = LHSmat;

    sel_nets = NaN(size(sel_nets_base,1),20); % format for probiotic layout
    sel_nets(:,probiotic_param_indexes) = repmat(null_probiotic,size(sel_nets,1),1);
    base_species_indexes = 1:20;
    base_species_indexes(probiotic_param_indexes) = [];
    sel_nets(:,base_species_indexes) = sel_nets_base;

    % In silico patient population
    sel_idx = 1:size(sel_nets,1); % when input workspace is same as used workspace, this is just a list of 1:number Virtual Patients
    POPinfo.virtualPatientParameters = sel_nets;
    POPinfo.virtualPatientComposition = StbleSS(sel_idx);
    
    % DEFAULT PROBIOTIC SIMULATION SETTINGS
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Set up run information to be called in simulate_CST_probiotic_object.m by 
    % defining the simulation components (POPinfo, PROBinfo, SIMinfo)
    
    % ~~~~~~~~~~ Acute Probiotic Dosing ~~~~~~~~~~ 
    PROBinfo.dose = 7.5; % Probiotic dose
    PROBinfo.numberDoses = 7; % Total doses given
    PROBinfo.doseFrequency = 1; % Every X days
    

    % ~~~~~~~~~~ Format Inputs for Different Simulation Types ~~~~~~~~~~ 
    if simulationType == "base" 
        % Acute Probiotic Dosing (no combination w/ antibiotic)
        SIMinfo.startProbiotic = 18; % Probiotic start point
        SIMinfo.stopSimulation = 365; % Simulation end point post probiotic treatment
        SIMinfo.initialComposition = displace_SS; % initial composition
        SIMinfo.alteredParamIndex = probiotic_param_indexes; % empty if no strain level changes
        SIMinfo.paramAlteration = null_probiotic; % null probiotic (average values)
        SIMinfo.paramAlterationType = 'absolute'; % ('None','fold addition', 'absolute')
        SIMinfo.ODEoptions = odeOptions;
        SIMinfo.ODEsolver = odeSolver;
        
        ABXinfo = NaN; % does not apply
    elseif simulationType == "combo" 
        % Acute Probiotic Dosing (no combination w or w/o antibiotic)
        PROBinfo.strainCharIndexes =  probiotic_param_indexes; % probiotic parameter indices
        PROBinfo.strainCharacteristics = null_probiotic; % null probiotic (no interaction terms)
       
        POPinfo.virtualPatientParameters = sel_nets; % in silico patients
        POPinfo.Indexes = sel_idx; % reference index for in silico paptients
        POPinfo.virtualPatientComposition = StbleSS(sel_idx);
        POPinfo.initialComposition = displace_SS; % displace from equilibrium
        POPinfo.odeSolver = odeSolver;
        POPinfo.odeOptions = odeOptions;

        ABXinfo = []; % define within run script
        SIMinfo = NaN; % does not apply
    elseif simulationType == "comboLV"
        % Custom Probiotic Dosing (w/ antibiotic, Cohen et a. 2020)
        
        % DEFINE DOSING REGIIMEN
        % Lactin-V was dosed as follows:
        %   Week 0: 5 days intravaginal MTZ
        %   Week 1: 4 days probiotic
        %   Weeks 2-11: 2/wk probiotic
    
        ABXinfo.addABX = true; % true: administer ABX, false: no ABX
        ABXinfo.dose = -3.82; % decay rate
        ABXinfo.duration = 5; % days
        
        % probiotic dosing regimen
        wk1 = [1 1 1 1 0 0 0]; % wk 1: 1/day 4 days
        wks = 10; % 2/wk for 10 wks
        int_lbs = repmat([0 1 0 0 1 0 0],1,wks);
        customDosing = find([wk1,int_lbs]); % get numerical values
    
        PROBinfo.strainCharIndexes =  probiotic_param_indexes; % probiotic parameter indices
        PROBinfo.strainCharacteristics = null_probiotic; % null probiotic (no interaction terms)
        PROBinfo.custom = true;
        PROBinfo.customProbDosing = customDosing;
        PROBinfo.numberDoses = 1; % days
        PROBinfo.doseFrequency = NaN; % every XX days
        PROBinfo.customEvalPoints = [1:7:84, 77+30, 77+90, 77+180];
    
        % NOTE: the virtualPatientParameters and Composition will be different
        % for the LV population
        POPinfo.initialComposition = displace_SS; % displace from equilibrium
        POPinfo.odeSolver = odeSolver;
        POPinfo.odeOptions = odeOptions;
        POPinfo.Indexes = sel_idx;

        SIMinfo = NaN; % does not apply
    else
        warning('simulationType must be base, combo, or comboLV')
    end
end