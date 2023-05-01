%% [MultiStable_NO,MonoStable_NO,CSTlabel,SSnum,rSSnum] = getNO_community_info(LHSmat,icArray,tspan,options)
%
%
%
function [MultiStable_NO,MonoStable_NO,CSTlabel,SSnum,rSSnum] = getNO_community_info(LHSmat,StbleSS,icArray,tspan,options)
    N = size(LHSmat,1);
    num_sp = 4;

    SSnum = NaN(size(icArray,1),N,num_sp);
    rSSnum = NaN(size(icArray,1),N,num_sp);

    CSTlabel = NaN(size(icArray,1),N);
    SStrack = NaN(size(icArray,1),N);
    SSnms = {'[nAB] CST-IV','[Li] CST-III','[oLB] CST-I/II/V','NoStbSS'};

    SSref = [0.911801263	0.059240659	0.028958078;
        0.146456352	0.758938299	0.09460535;
        0.153186234	0.09522381	0.751589956];

    deltT = 10;

    for ic_id = 1:size(icArray,1)
        for set_id = 1:N
            params = LHSmat(set_id,:);

            scle = StbleSS{set_id}{1};
            if ~isempty(scle)
                y0 = icArray(ic_id,:)*sum(scle(1,:));
                [t, y] = ode45(@lhs_ode_gLV,tspan,y0,options,params);
                yend = y(end,:);
            end

            if t(end) < tspan(end) - deltT
                yend = NaN(1,num_sp);
            else
                stmp = sum((y(end,:) - y(size(y,1) - deltT,:)).^2);
                SStrack(ic_id,set_id) = stmp;
                disp(stmp > 1)
            end

            SSnum(ic_id,set_id,:) = yend;

            tmp = yend/sum(yend);
            rSSnum(ic_id,set_id,:) = tmp;

            % Get CST Label
            if sum(isnan(tmp)) == 0
                dif = sum((SSref - repmat(tmp(1:3),size(SSref,1),1)).^2,2);
                [val,id] = min(dif);
                CSTlabel(ic_id,set_id) = id;
            else
                CSTlabel(ic_id,set_id) = 4;
            end

            disp(strcat("IC: ", num2str(icArray(ic_id,:)), "Set #:", num2str(set_id)))
        end
    end

    %% SpotCheck Sets

    set_id = 2;
    squeeze(SSnum(:,set_id,:))

    %% Number of NO sets:
    NO_CST = CSTlabel == 1;
    num_NO_CST = sum(sum(NO_CST) > 0);
    disp(strcat("There were ", num2str(num_NO_CST), " sets that are NO w/o probiotic"))

    CST_sel_id = sum(NO_CST) > 0;

    numberSS = NaN(N,1);
    for ms_id = 1:N
        tmp = CSTlabel(:,ms_id);
        numUniq = length(unique(tmp));
        numberSS(ms_id) = numUniq;
    end

    multi_stable_NO_num = sum(numberSS' > 1 & sum(NO_CST) > 0);

    disp(strcat("There were ", num2str(multi_stable_NO_num), " sets that were multi-stable"))

    %% NO parameters idx
    MultiStable_NO = numberSS' > 1 & sum(NO_CST) > 0;
    MonoStable_NO = numberSS' == 1 & sum(NO_CST) > 0;

end


    