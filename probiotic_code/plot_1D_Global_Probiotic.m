function plot_1D_Global_Probiotic(fdr_nm)

%%
    [param_names,SSnms,sp_cols,~] = get_naming_terms();
    lst = dir(fdr_nm);
    ws_nms = string({lst.name});
    
    params_to_view = [8,12,16,17,18,19];

    num_iters = sum(contains(ws_nms,'.mat'))/length(params_to_view);
    
    NoResponse = NaN(length(params_to_view),num_iters);
    LiResponse = NaN(length(params_to_view),num_iters);
    oLBResponse = NaN(length(params_to_view),num_iters);
    PResponse = NaN(length(params_to_view),num_iters);
    pnm = extractBefore(ws_nms,'_');
    %%
    for param_id = 1:length(params_to_view)
        
        idx_ws = pnm == string(params_to_view(param_id));
        sel_ws = ws_nms(idx_ws);
        selpm = param_names(params_to_view(param_id));
        all_count = NaN(4,length(sel_ws));
        all_perc = NaN(4,length(sel_ws));
        titles = string(NaN(length(sel_ws),1));
        magorder = NaN(1,length(sel_ws));
        
        % Pull data from workspaces
        for i = 1:length(sel_ws)
            load(strcat(fdr_nm,sel_ws(i)),'SSnms','CST_post','warnall','pmag')
            CST_post_rmerr = CST_post(warnall == "");
            TabPost = tabulate(SSnms(CST_post_rmerr));
            disp(strcat(" For: x + ", num2str(pmag),"x"))
            disp(TabPost)
            count_out = NaN(4,1);
            perc_out = NaN(4,1);
            for k = 1:length(SSnms)
                id = sum(CST_post_rmerr == k);
                count_out(k) = id;
                perc_out(k) = id/length(CST_post_rmerr);
            end
            all_count(:,i) = count_out;
            all_perc(:,i) = perc_out;
            titles(i) = num2str(pmag);
            magorder(i) = pmag;
        end
        
        % Prepare data for plotting
        [magvals,id_sort] = sort(magorder);
        dat_in = all_perc;
        
        Y = dat_in;
        Y = bsxfun(@rdivide, Y, sum(Y,1));
        
        subplot(2,3,param_id)
        area(magvals,Y(:,id_sort)')
%         xlim([-4,4])
        colororder(sp_cols)
    
        xlabel('Fold Addition')
        ylabel('Percent Samples')
        title(param_names(params_to_view(param_id)))
        legend(SSnms)
        ylim([0 1])
        
    
    %   Label with specific values picked
    %     ax = gca;
    %     xticklabels(magvals)
    %     ax.XTick = magvals;
    
        NoResponse(param_id,:) = dat_in(1,id_sort);
        LiResponse(param_id,:) = dat_in(2,id_sort);
        oLBResponse(param_id,:) = dat_in(3,id_sort);
        PResponse(param_id,:) = dat_in(4,id_sort);
        set(gca,'fontsize',14)
    end
end