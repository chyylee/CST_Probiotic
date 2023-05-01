%% Fig4_Process_Systematic_LactinV_Result.m
%
% GOAL: 
%  1. Understand what parameters could potentially explain the results
%       observed in the Lactin-V study (Box-plots, Volcanos)
%  2. Understand what parameters would make a probiotic more effective for
%       a Lactin-V like regimen (PLS-DA)
%           - Requires PLS Toolbox and MATLAB 2017
%
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% 1. Load workspace generated in Fig4_Systematic_LactinV_Probiotic.m
output_fdr_nm = 'result_workspaces/F4/';
load(strcat(output_fdr_nm,'F4_Systematic_ShortTerm_Probiotic.mat'))
    % all_select_outcomes: Numbers Strains x Simulated Patients x Time
    %       Points x [time point, nAB, Li, oLB, Probiotic] Abundance
    % all_warnall: Indicates runs without errors (1 = no error, 0 = error)
    % strainListLHS: Probiotic strain parameter values (Numbers Strains x
    %       Strain Parameters)
    % POPinfo: Structure with information on simulated population
    % PROBinfo: Structure with information about probiotic
    % SIMinfo: Structure with information about simulation (length, ODE
    %       solver conditions, etc.)


%% 2. Pull Data of Strains at Lactin-V Evaluation Points
% NOTE: Lactin-V is evaluated at 12 weeks and 24 weeks after initial
% probiotic dosing. This corresponds to index 12 and 14 of
% all_select_outcomes

[param_names,SSnms,sp_cols,time_names] = get_naming_terms();
time_id = [5 8]; % 1mo or 12mo
numPatients = size(all_select_outcomes,2);
numStrains = size(strainListLHS,1);
numODEwarnings = sum(all_warnall==0,2);

errorThresh = 1; % must have a least one successful run (no ODE solver errors)

% collect information about response counts and frequencies
strainCounts1mo = NaN(numStrains,4);
strainFrequency1mo = NaN(numStrains,4);
strainCounts12mo = NaN(numStrains,4);
strainFrequency12mo = NaN(numStrains,4);

% loop through each "strain" (strainListLHS combination)
for strain_id = 1:numStrains
    warn_id = all_warnall(strain_id,:) == 1;
    if sum(warn_id) > errorThresh
        CST_post = get_CST_post(squeeze(all_select_outcomes(strain_id,warn_id,time_id(1),2:end)));
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts);
        strainCounts12(strain_id,:) = Counts;
        strainFrequency1mo(strain_id,:) = Frequency;

        CST_post = get_CST_post(squeeze(all_select_outcomes(strain_id,warn_id,time_id(2),2:end)));
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts);
        strainCounts12mo(strain_id,:) = Counts;
        strainFrequency12mo(strain_id,:) = Frequency;
    end
end
disp('Maximum Response Frequency 1mo: nAB, Li, oLB, Probiotic:')
disp(max(strainFrequency1mo))

disp('Maximum Response Frequency 12mo: nAB, Li, oLB, Probiotic:')
disp(max(strainFrequency12mo))


errorThresh = 500; % threshold for the number of successful runs (244 gives 1000 "good" runs)


%% 6. RUN PLS-DA
X = strainFrequency1mo; % or strainFrequency24
[~,lb] = max(X > 0.5,[],2);
clus_ssize = sum(dummyvar(lb));
lb(isnan(X(:,1))) = -1;

prism_dat = [];
Xblock = [];
yblock = [];
c = 1;
for i = [1 2 3 4]
    tmp = strainListLHS(lb == i & numODEwarnings <= errorThresh,:);
    Xblock = [Xblock;tmp];
    yblock = [yblock;repmat(c,size(tmp,1),1)];
    c = c + 1;
end

yblock(yblock == 3) = 2;
yblock(yblock == 4) = 2;
classes = {'nAB','LB'};

% classes = {'nAB','Li','oLB','P'};
ynames = param_names(PROBinfo.strainCharIndexes);
Yblock = dummyvar(yblock);
filettl = 'XXXXXXXXXXX-1mo-ST-BVefficacy';
ttl = '';
yname = '';
ep = 1;
[plsda_model,cls_error,cls_error_CV,ellip_points] = fPLSDA(Xblock,...
    Yblock,ynames,filettl,classes,yname,ttl,1,ep);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% For copying andd pasting into PRISM (XY DATA)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X = NaN(size(LV1_scores,1),length(classes)+1); % compile scores by response type
X(:,1) = LV1_scores;
E = NaN(size(ellip_points,2)*length(classes),length(classes)+1); % points to draw ellipses
c = 1;
for i = 1:length(classes)
    X(yblock(:,i) == 1,i+1) = LV2_scores(yblock(:,i) == 1);
    ep = c + size(ellip_points,2)-1;

    E(c:ep,1) = squeeze(ellip_points(i,:,1));
    E(c:ep,i+1) = squeeze(ellip_points(i,:,2));
    c = c + size(ellip_points,2);

end

% Re-order the PLS Loadings
[sorted_LVloads,i] = sort(LV2_loads);

%%
function plot_pca(score,coeff,explained,pcn1,pcn2,xnames,Y,classes)

    subplot(1,2,1)
    if isempty(classes)
        p = scatter(score(:,pcn1), score(:,pcn2),40,Y,'filled');
        colorbar()
        hl = refline([0,0]);
        hl.Color = 'k';
        vl = xline(0);
    else
        p = gscatter(score(:,pcn1), score(:,pcn2), Y);
        
        hl = refline([0,0]);
        hl.Color = 'k';
        vl = xline(0);
        legend([p],classes)
    end

    xlabel(['PC 1 - ', num2str(round(explained(pcn1),2)), '%'])
    ylabel(['PC 2 - ', num2str(round(explained(pcn2),2)), '%'])
    set(gca,'fontsize',14)
    title('Scores Plot')


%     subplot(1,2,2)
%     scatter(coeff(:,pcn1),coeff(:,pcn2))
%     text(coeff(:,pcn1),coeff(:,pcn2), xnames)
%     xlabel(['PC 1 - ', num2str(round(explained(pcn1),2)), '%'])
%     ylabel(['PC 2 - ', num2str(round(explained(pcn2),2)), '%'])
%     hl = refline([0,0]);
%     hl.Color = 'k';
%     vl = xline(0);
%     set(gca,'fontsize',14)
%     title('Coefficient Plot')
end
