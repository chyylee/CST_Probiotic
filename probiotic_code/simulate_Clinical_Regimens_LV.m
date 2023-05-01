%% [select_outcomes,  meanstdDuring, warnall, nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo)
%
% GOAL:
%
% Assumes in silico population is BV+
%
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function [select_outcomes,  meanstdDuring, warnall, nonBVflag] = simulate_Clinical_Regimens_LV(ABXinfo,PROBinfo,POPinfo)
    inSilicoPopSize = size(POPinfo.virtualPatientParameters,1);
    warnall = string(NaN(inSilicoPopSize,1));
    nonBVflag = NaN(inSilicoPopSize,1);
    
    eq_time = 5;
    st_pr = eq_time + ABXinfo.duration;
    
    if PROBinfo.custom 
        dose_lbls = PROBinfo.customProbDosing;
        stop_point = dose_lbls(end); 
        extraEvalPoints = st_pr + PROBinfo.customEvalPoints;
    elseif PROBinfo.Number > 0
        dose_lbls = st_pr:PROBinfo.doseFrequency:st_pr+PROBinfo.doseFrequency*(PROBinfo.numberDoses-1);
        stop_point = dose_lbls(end); 
        extraEvalPoints = [];
    else
        stop_point = st_pr; 
        extraEvalPoints = [];
    end

    eval_points = [extraEvalPoints,...
        (stop_point - st_pr)/2, ...
        stop_point, stop_point + 7, ...
        stop_point + 14, ...
        stop_point + 28, ...
        stop_point + 12*7, ...
        stop_point + 24*7, ...
        stop_point + 52*7];
    
    select_outcomes = NaN(inSilicoPopSize,length(eval_points),5);
    meanstdDuring = NaN(inSilicoPopSize,2,4); % mean and stand deviation during treatment
  
    parfor net_id = 1:inSilicoPopSize
        disp(strcat("RUNNING ... #", num2str(net_id)))
        obsSS = POPinfo.virtualPatientComposition{net_id}{:};% use predicted absolute abundance for community
        [~,cd] = max(obsSS(:,1)); % finds BV+ state if there are multiple
        ybase = obsSS(cd,:); 
        
        yscale = POPinfo.initialComposition;
        y0 = [ybase,0] + sum(ybase).*yscale;
        
        % parameters values
        params = POPinfo.virtualPatientParameters(net_id,:);
        params(PROBinfo.strainCharIndexes) = PROBinfo.strainCharacteristics;
        prob_amt = PROBinfo.dose;
        
        tspanS = [0 eq_time];
        % run base simulation
        [ts, ys] = POPinfo.odeSolver(@lhs_ode_gLV,tspanS,y0,POPinfo.odeOptions,params);
        [warnmsg1, ~] = lastwarn; % collect any ode solver warnings pre-probiotic
        
        checkInitialState = ys(end,1)./sum(ys(end,:));
        if checkInitialState < 0.8
            nonBVflag(net_id) = 1;
        end
        
        % run the ABX simulation
        
        if ABXinfo.addABX
            tspanA = [eq_time st_pr];
            paramsABX = params;
            paramsABX(1) = params(1) + ABXinfo.dose; % assumes growth rate of nAB is the first parameter
            [tabx, yabx] = POPinfo.odeSolver(@lhs_ode_gLV,tspanA,ys(end,:),POPinfo.odeOptions,paramsABX);
            yprob = yabx(end,:);
        else
            tabx = [];
            yabx = [];
            yprob = ys(end,:);
        end
            
            
        % spike in probiotic
        tcol = [ts;tabx];
        ycol = [ys;yabx];
        tlast = eval_points(end);

            
        if PROBinfo.numberDoses > 0
            yprob(4) = PROBinfo.dose;
            
            dose_sep = dose_lbls(2:end) - dose_lbls(1:end-1); % get dose spacing
            
            st_counter = tcol(end);
            % Repeat dosing if needed
            for prob_iter = 1:length(dose_lbls)-1
                tspanP = [st_counter, st_counter+dose_sep(prob_iter)];
                st_counter = st_counter+dose_sep(prob_iter);
            
                [tp, yp] = POPinfo.odeSolver(@lhs_ode_gLV,tspanP,yprob,POPinfo.odeOptions,params);
            
                tcol = [tcol; tp];
                ycol = [ycol; yp];
            
                % Set up initial condition for next dose
                yprob = yp(end,:);
                yprob(4) = yp(end,4)+prob_amt;
            end
            tspanP = [st_counter,st_counter + tlast];
            [tp, yp] = POPinfo.odeSolver(@lhs_ode_gLV,tspanP,yprob,POPinfo.odeOptions,params);
            
            tcol = [tcol; tp];
            ycol = [ycol; yp];
        else
           tspanP = [st_pr, st_pr + tlast];
           [tp, yp] = POPinfo.odeSolver(@lhs_ode_gLV,tspanP,yprob,POPinfo.odeOptions,params);
           tcol = [tcol; tp];
           ycol = [ycol; yp];

        end
  
            
        % collect any ode solver warnings post-probiotic
        [warnmsg2, ~] = lastwarn; 

        % save abundances at evaluation points
        select_outcomes(net_id,:,:) = pull_time_series_data(eval_points,tcol,ycol);

        % save profile during probiotic dosing
        if PROBinfo.dose > 0
            [~,idx_probiotic_start] = min(abs(tcol - st_pr)); % finds closest evaluation time point
            [~,idx_probiotic_stop] = min(abs(tcol - stop_point)); % finds closest evaluation time point
    
            meanstdDuring(net_id,:,:) = [mean(ycol(idx_probiotic_start:idx_probiotic_stop,:));
                std(ycol(idx_probiotic_start:idx_probiotic_stop,:))];
        else
            meanstdDuring = NaN;
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