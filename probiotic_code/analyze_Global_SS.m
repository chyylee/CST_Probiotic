%% analyze_Global_SS.m
% analyze_Global_SS(ws_nm)
%
% Can be used to assess the steady state configurations for the output of
% the global sensitivity analysis LHS parameter sets.
%
% INPUT:
%   * ws_nm: name of LHS workspace generated by the global sensitivity
%       analysis
%
% OUTPUT:
%   * New workspace with the saved steady-state configurations and a plot
%   of the frequency of each.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Christina Y. Lee
% University of Michigan
% Jan 12, 2022
% Update: Jan 20, 2023 (removed loading analytic mat, added symbolic_solns
% call)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function analyze_Global_SS(LHSmat,sp_names)

    % Generate parameter & model information
    N = length(sp_names); % Number of species in model
    [S, Jmat] = symbolic_solns(N); % symbolic solutions to N species model
    
    % Parameter names given input species names
    [nm_out1] = generate_parameter_names(sp_names);
    [nm_out2] = generate_coeff_labels('\alpha',sp_names);
    param_names = horzcat(nm_out1{1:length(sp_names)},nm_out2);


    %% Caculate steady-states for a given parameter set in LHSmat
    StbleSS = cell(size(LHSmat,1),1);
    ALLSS = cell(size(LHSmat,1),1);
    eigSS = cell(size(LHSmat,1),1);

    parfor i = 1:size(LHSmat,1) % loops through "simulated samples"
        disp(strcat("Set #", num2str(i)))
        params = LHSmat(i,:);
        [StableStates,Sval,eigvals,~] = calc_SS_stability(N,params,S,Jmat);
        StbleSS{i} = {StableStates};
        ALLSS{i} = {Sval};
        eigSS{i} = {eigvals};
    end


    %% 3) ANALYZE MODEL SS MONO- VS MULTI- STTABILITY

    % Names for possible steady-states in model
    poss_SSnames = {'Co-elim','NO','Li','oLB',...
            'NO&Li', 'NO&oLB', 'Li&oLB', 'NOLi&oLB'};
    poss_SS = strings(size(StbleSS,1),1);

    % Loop tthrough and name each state (looks to see how many states are
    % possible, 0: Unstable, 1: Mono-stable, 2: Bi-stable... etc.)
    for i = 1:size(StbleSS,1)
        temp = StbleSS{i}{:};
        num_SS = size(temp,1);
        if num_SS == 0
            poss_SS(i) = "All Unstable";
        else
            flg = strcat(num2str(num_SS),'SS');
            fnm = [];
            for j = 1:num_SS
                stmp = temp(j,:) > 0;
                nm = strcat(sp_names{stmp});
                if sum(stmp) > 1
                    nm = strcat(nm," or ");
                elseif sum(stmp) == 1
                    nm = strcat(nm," or ");
                else
                    nm = strcat('',nm);
                end
                fnm = strcat(fnm,nm);
            end
            cfnm = char(strcat(flg,": ",fnm));
            ffnm = cfnm(1:length(cfnm) - 4);
            poss_SS(i) = string(ffnm);
        end
    end

    tabulate(poss_SS) % countts response pattern (prior to removing unstable)

    %% 4) COMPILE OUTPUT
    noUnstable = poss_SS; % remove unstable states
    noUnstable(contains(poss_SS,"Unstable")) = [];

    numUS = length(poss_SS) - length(noUnstable);

    out = tabulate(noUnstable);
    SS_names = out(:,1);
    SS_counts = [out{:,2}]';
    SS_percent = [out{:,3}]';

    idx1SS = contains(SS_names,'1SS');

    monosum = sum(SS_percent(idx1SS));
    multisum = sum(SS_percent(~idx1SS));

    % PLOT REULT
    [v,i] = sort(SS_percent,'descend');
    n = length(v);

    figure(1)
    bar(v)
    xticks(1:n)
    xticklabels(SS_names(i))
    xtickangle(270)
    hold on
    text(1:n,v,string(round(v,3)'),'vert','bottom','horiz','center'); 
    title(strcat(" SS-Configs (N = ", num2str(sum(SS_counts)),")",...
        " (",num2str(numUS), " no stable state)"))
    ylabel('Percent of Stable Steady-States')

    dim = [.65 .8 .2 .1];
    str = strcat("Mono-Stable: ", num2str(round(monosum,2))...
        , "%") + newline + strcat("Multi-Stable: ", num2str(round(multisum,2)), "%");
    annotation('textbox',dim,'String',str)
    
    colors = [255 242 204;
        193 147 171;
        108 67 197;
        180 198 231;
        139 23 24;
        194 194 194;
        21 225 204;
        168 99 255;
        174 164 255;
        255 183 255;
        255 244 111;
        231 0 255;
        0 2 255;
        224 170 150;
        80 15 201;
        114 157 195;
        169 208 142;
        84 130 53;
        55 86 35;
        191 143 0;
        0 0 0]./255;

    %% Save File
  out_nm = strcat('SSConfig-Analysis.mat');

  save(out_nm,'StbleSS','S','Jmat','N','poss_SS','noUnstable',...
        'SS_names','SS_counts','SS_percent','monosum','multisum','numUS','poss_SSnames',...
        'LHSmat','colors','sp_names','param_names','ALLSS','eigSS')
  if N == 3
      figure(2)
      analyze_Global_CST_SS(out_nm)
  end
end

function analyze_Global_CST_SS(fn)
    load(fn)
    all_nm_CST = [];
    for idx = 1:length(StbleSS)
        run_mat = cell2mat(StbleSS{idx});
        [nmf] = get_VALENCIA_class(run_mat); % CONVERT TO VALENCIA CST TYPE
        all_nm_CST = [all_nm_CST; nmf];
    end

    tabulate(all_nm_CST)

    % PLOT THE OUTPUT WITH NO-SS REMOVED
    poss_SSv = all_nm_CST;

    noUnstablev = poss_SSv; % remove unstable states
    noUnstablev(contains(poss_SSv,"0SS")) = [];

    numUS = length(poss_SSv) - length(noUnstablev);

    out = tabulate(noUnstablev);
    SS_names_CST = out(:,1);
    SS_counts_CST = [out{:,2}]';
    SS_percent_CST = [out{:,3}]';

    idx1SS = contains(SS_names_CST,'1SS');

    monosum = sum(SS_percent_CST(idx1SS));
    multisum = sum(SS_percent_CST(~idx1SS));

    [v,i] = sort(SS_percent_CST,'descend');
    n = length(v);

    bar(v)
    xticks(1:length(v))
    xticklabels(SS_names_CST(i))
    xtickangle(270)
    hold on

    text(1:n,v,string(round(v,3)'),'vert','bottom','horiz','center'); 
    ylabel('Percent of Stable Steady-States')

    dim = [.65 .8 .2 .1];
    str = strcat("Mono-Stable: ", num2str(round(monosum,2))...
        , "%") + newline + strcat("Multi-Stable: ", num2str(round(multisum,2)), "%");
    annotation('textbox',dim,'String',str)
    title('Equilibrium Behaviors (Model SS Classified to CSTs)')


    % 4. Append VALCENIA CST to data
    save(fn,'all_nm_CST','SS_names_CST','SS_percent_CST','SS_counts_CST',...
        '-append')
end