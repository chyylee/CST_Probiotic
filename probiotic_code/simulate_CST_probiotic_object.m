%% [select_outcomes,nonBVflag,tall,yall,warnall,dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo)
%
% Runs probiotic simulations for the following scenarios:
%   - A given parameter set (probiotic included)
%   - Modification to each parameter set
%       - fold addition
%       - absolute change (set to new value)
%
% --------------------- EXAMPLE RUNS START --------------------- 
% Note: sel_nets, sel_idx, and StbleSSPre will need to be loaded from the
% global sensitivity analysis that generated the virtual patient
% populations
%
% Example run with no parameter change
% POPinfo.virtualPatientParameters = sel_nets;
% POPinfo.virtualPatientComposition = StbleSSPre(sel_idx);
% 
% PROBinfo.dose = 7.5;
% PROBinfo.numberDoses = 7;
% PROBinfo.doseFrequency = 1;
% 
% SIMinfo.startProbiotic = 5;
% SIMinfo.stopSimulation = 30;
% SIMinfo.initialComposition = [0.98 0.01 0.01 0];
% SIMinfo.alteredParamIndex = []; % empty if no strain level changes
% SIMinfo.paramAlterationType = 'None'; % ('None','fold addition', 'absolute')
% SIMinfo.paramAlteration = [];
% SIMinfo.ODEoptions = [];
% SIMinfo.ODEsolver = @ode15s;
% 
% [select_outcomes,nonBVflag,tall,yall,warnall, ...
%     dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo)
% % 
% % Check with perturbation
% POPinfo.virtualPatientParameters = sel_nets;
% POPinfo.virtualPatientComposition = StbleSSPre(sel_idx);
% 
% PROBinfo.dose = 7.5;
% PROBinfo.numberDoses = 7;
% PROBinfo.doseFrequency = 1;
% 
% SIMinfo.startProbiotic = 5;
% SIMinfo.stopSimulation = 30;
% SIMinfo.initialComposition = [0.98 0.01 0.01 0];
% SIMinfo.alteredParamIndex = [1 2 3]; % empty if no strain level changes
% SIMinfo.paramAlterationType = 'absolute'; % ('None','fold addition', 'absolute')
% SIMinfo.paramAlteration = [-0.5 0 -.5];
% SIMinfo.ODEoptions = [];
% SIMinfo.ODEsolver = @ode15s;
% 
% [select_outcomes,nonBVflag,tall,yall,warnall,...
%       dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo);
% --------------------- EXAMPLE RUNS END --------------------- 
%%
function [select_outcomes,nonBVflag,tall,yall,warnall,dose_lbls] = simulate_CST_probiotic_object(POPinfo,PROBinfo,SIMinfo)
   %%
    NumPatients = size(POPinfo.virtualPatientParameters,1);
    NumSpecies = length(SIMinfo.initialComposition);

    warnall = string(NaN(NumPatients,1));
    tall = cell(NumPatients,1);
    yall = cell(NumPatients,1);
    nonBVflag = NaN(NumPatients,1);

    st_pr = SIMinfo.startProbiotic;
    dose_sep = PROBinfo.doseFrequency;
    dose_rep = PROBinfo.numberDoses;
    dose_lbls = st_pr:dose_sep:st_pr+dose_sep*(dose_rep-1);
    
    stop_point = dose_lbls(end); % last dose
    eval_points = [(stop_point - st_pr)/2, ...
        stop_point, stop_point + 7, ...
        stop_point + 14, ...
        stop_point + 28, ...
        stop_point + 12*7, ...
        stop_point + 24*7, ...
        stop_point + 52*7];

   select_outcomes = NaN(NumPatients,length(eval_points),NumSpecies+1);

   StbleSSPre = POPinfo.virtualPatientComposition;
   sel_nets = POPinfo.virtualPatientParameters;
   yscale = SIMinfo.initialComposition;
   simulationType = SIMinfo.paramAlterationType;
   pidx = SIMinfo.alteredParamIndex;
   pmag = SIMinfo.paramAlteration;
   prob_dose = PROBinfo.dose;
   odeSolver = SIMinfo.ODEsolver;
   options = SIMinfo.ODEoptions;
   tlast = SIMinfo.stopSimulation;
   
   possibleSimulationTypes = ["None","Fold Addition","Absolute"];
   codedSimulationType = lower(string(simulationType)) == lower(possibleSimulationTypes);
   if sum(codedSimulationType) == 0
       disp('WARNING: Please select None, Fold Addition, or Absolute')
   end
   codedSimulationType = find(codedSimulationType);
   parfor net_id = 1:NumPatients
        % Set up initial conditions and parameters
        obsSS = StbleSSPre{net_id}{:};% use predicted absolute abundance for community
        [~,cd] = max(obsSS(:,1)); % get state with max BV bacteria
        ybase = obsSS(cd,:); 
%         y0 = sum(ybase).*yscale;
        y0 = [ybase,0] + sum(ybase).*yscale; % update 04/19/2023

        % parameters values
        base_params = sel_nets(net_id,:);
        
        % here add option about parameter change
        if codedSimulationType == 1
            params = base_params;
        elseif codedSimulationType == 2
            params = base_params;
            for k = 1:length(pidx)
                params(pidx(k)) = base_params(pidx(k)) + abs(base_params(pidx(k)))*pmag(k);
            end
        elseif codedSimulationType == 3
            params = base_params;
            params(pidx) = pmag;
        else
            disp('Warning: incorrect input for simulation type')
        end
        
        prob_amt = prob_dose;
        
        tspanS = [0 st_pr];
        % run base simulation
        [ts, ys] = odeSolver(@lhs_ode_gLV,tspanS,y0,options,params);
        [warnmsg1, ~] = lastwarn; % collect any ode solver warnings pre-probiotic
        
        rys = ys(end,1)./sum(ys(end,:));
        if rys < 0.8
            nonBVflag(net_id) = 1;
        end

        % spike in probiotic
        yprob = ys(end,:);
        yprob(4) = prob_amt;

        tcol = ts;
        ycol = ys;
        % Repeat dosing if needed
        for prob_iter = 1:dose_rep
            if prob_iter == 1 && prob_iter == dose_rep
                tspanP = [st_pr, st_pr + tlast];
            elseif prob_iter == 1 && prob_iter ~= dose_rep
                tspanP = [st_pr, st_pr + dose_sep];
                nstr_pr = st_pr + dose_sep;
            elseif prob_iter == dose_rep
                tspanP = [nstr_pr, nstr_pr + tlast];
            else
                tspanP = [nstr_pr, nstr_pr + dose_sep];
                nstr_pr = nstr_pr + dose_sep;
            end
            
            [tp, yp] = odeSolver(@lhs_ode_gLV,tspanP,yprob,options,params);

            tcol = [tcol; tp];
            ycol = [ycol; yp];

            % Set up initial condition for next dose
            yprob = yp(end,:);
            yprob(4) = yp(end,4)+prob_amt;
        end

        tmp = NaN(8,5);
        for i = 1:length(eval_points)
            [~,idx] = min(abs(tcol - eval_points(i)));
            if ~isempty(idx)
                tmp(i,:) = [tcol(idx(1)),ycol(idx(1),:)];
            else
                tmp(i,:) = NaN(1,5);
            end
        end
        select_outcomes(net_id,:,:) = tmp;

        % collect any ode solver warnings post-probiotic
        [warnmsg2, ~] = lastwarn; 

        % save ODE run information
        tall(net_id) = {[tcol]};
        yall(net_id) = {[ycol]};

        % keep track of runs with and without warningss
        warnall(net_id) = string(strcat(warnmsg1,warnmsg2)); 
        lastwarn('') % clear warnings for next parameter set
   end
end