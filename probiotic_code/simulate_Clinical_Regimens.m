function [select_outcomes,  meanstdDuring, warnall, nonBVflag] = simulate_Clinical_Regimens(ABXinfo,PROBinfo,POPinfo)
    inSilicoPopSize = size(POPinfo.virtualPatientParameters,1);
    warnall = string(NaN(inSilicoPopSize,1));
    nonBVflag = NaN(inSilicoPopSize,1);
    
    eq_time = 5;
    st_pr = eq_time + ABXinfo.duration;
    

    doseFrequency = PROBinfo.doseFrequency;
    numberDoses = PROBinfo.numberDoses;
    if PROBinfo.numberDoses ~= 0
        dose_lbls = st_pr:doseFrequency:st_pr+doseFrequency*(numberDoses-1);
    else
        dose_lbls = st_pr;
    end

    stop_point = dose_lbls(end); % last dose

    eval_points = [(stop_point - st_pr)/2, ...
        stop_point, stop_point + 7, ...
        stop_point + 14, ...
        stop_point + 28, ...
        stop_point + 12*7, ...
        stop_point + 24*7, ...
        stop_point + 52*7];
    
    select_outcomes = NaN(inSilicoPopSize,length(eval_points),5);
    meanstdDuring = NaN(inSilicoPopSize,2,4); % mean and stand deviation during treatment
    virtualPatientComposition = POPinfo.virtualPatientComposition;
    initialComposition = POPinfo.initialComposition;
    virtualPatientParameters = POPinfo.virtualPatientParameters;
    strainCharIndexes = PROBinfo.strainCharIndexes;
    strainCharacteristics = PROBinfo.strainCharacteristics;
    odeSolver = POPinfo.odeSolver;
    probioticDose = PROBinfo.dose;
    odeOptions = POPinfo.odeOptions;
    addABX = ABXinfo.addABX;
    ABXdose = ABXinfo.dose;
    parfor net_id = 1:inSilicoPopSize
        disp(strcat("RUNNING ... #", num2str(net_id)))
        obsSS = virtualPatientComposition{net_id}{:};% use predicted absolute abundance for community
        [~,cd] = max(obsSS(:,1)); % finds BV+ state if there are multiple
        ybase = obsSS(cd,:); 
        
        yscale = initialComposition;
        y0 = [ybase,0] + sum(ybase).*yscale;
        
        % parameters values
        params = virtualPatientParameters(net_id,:);
        params(strainCharIndexes) = strainCharacteristics;
        prob_amt = probioticDose;
        
        tspanS = [0 eq_time];
        % run base simulation
        [ts, ys] = odeSolver(@lhs_ode_gLV,tspanS,y0,odeOptions,params);
        [warnmsg1, ~] = lastwarn; % collect any ode solver warnings pre-probiotic
        
        checkInitialState = ys(end,1)./sum(ys(end,:));
        if checkInitialState < 0.8
            nonBVflag(net_id) = 1;
        end
        
        % run the ABX simulation
        if addABX
            tspanA = [eq_time st_pr];
            paramsABX = params;
            paramsABX(1) = params(1) + ABXdose; % assumes growth rate of nAB is the first parameter
            [tabx, yabx] = odeSolver(@lhs_ode_gLV,tspanA,ys(end,:),odeOptions,paramsABX);
            yprob = yabx(end,:);
        else
            tabx = [];
            yabx = [];
            yprob = ys(end,:);
        end
            
            
        % spike in probiotic
        yprob(4) = probioticDose;
        
        tcol = [ts;tabx];
        ycol = [ys;yabx];
        dose_rep = numberDoses;
        dose_sep = doseFrequency;
        tlast = eval_points(end);
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
        
            [tp, yp] = odeSolver(@lhs_ode_gLV,tspanP,yprob,odeOptions,params);
        
            tcol = [tcol; tp];
            ycol = [ycol; yp];
        
            % Set up initial condition for next dose
            yprob = yp(end,:);
            yprob(4) = yp(end,4)+prob_amt;
        end
  
            
        % collect any ode solver warnings post-probiotic
        [warnmsg2, ~] = lastwarn; 

        % save abundances at evaluation points
        select_outcomes(net_id,:,:) = pull_time_series_data(eval_points,tcol,ycol);

        % save profile during probiotic dosing
        [~,idx_probiotic_start] = min(abs(tcol - st_pr)); % finds closest evaluation time point
        [~,idx_probiotic_stop] = min(abs(tcol - stop_point)); % finds closest evaluation time point
        
        if prob_amt > 0
            meanstdDuring(net_id,:,:) = [mean(ycol(idx_probiotic_start:idx_probiotic_stop,:));
                std(ycol(idx_probiotic_start:idx_probiotic_stop,:))];
        else
           meanstdDuring(net_id,:,:) = NaN;
        end
        % keep track of runs with and without warningss
        warnall(net_id) = string(strcat(warnmsg1,warnmsg2)); 
        lastwarn('') % clear warnings for next parameter set
    end
end

function [select_out] = pull_time_series_data(eval_points,tcol,ycol)
    for i = 1:length(eval_points)
        [~,idx] = min(abs(tcol - eval_points(i))); % finds closest evaluation time point
        if ~isempty(idx)
            select_out(i,:) = [tcol(idx(1)),ycol(idx(1),:)];
        end
    end
end
