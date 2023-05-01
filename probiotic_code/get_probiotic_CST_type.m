%%

% Modifying "analyze_Global_CST_SS.m to work with probiotic.
% Method: summing probiotic with oLB relative abundance


%%
%%
function get_probiotic_CST_type(StbleSS,fn)
    %load(fn)
    all_nm = [];
    for idx = 1:length(StbleSS)
        run_mat0 = cell2mat(StbleSS{idx});
        run_mat = [run_mat0(:,1:2), run_mat0(:,3) + run_mat0(:,4)]; % add oLB and P
        [nmf] = get_VALENCIA_class(run_mat);
        all_nm = [all_nm; nmf];
    end

    tabulate(all_nm)

    % PLOT THE OUTPUT WITH NO-SS REMOVED
    poss_SS_CST = all_nm;

    noUnstablev = poss_SS_CST; % remove unstable states
    noUnstablev(contains(poss_SS_CST,"0SS")) = [];

    numUS = length(poss_SS_CST) - length(noUnstablev);

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
    set(gca,'fontsize',16)

    %%
    % 4. Append VALCENIA CST to data
    save(fn,'all_nm','SS_names_CST','SS_percent_CST','SS_counts_CST',...
        '-append')
end