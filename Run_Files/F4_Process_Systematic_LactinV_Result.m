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
load(strcat(output_fdr_nm,'F4_Systematic_LactinV_Probiotic.mat'))
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
time_id = [12 14]; % 12 or 24 weeks
numPatients = size(all_select_outcomes,2);
numStrains = size(strainListLHS,1);
numODEwarnings = sum(all_warnall==0,2);

errorThresh = 1; % must have a least one successful run (no ODE solver errors)

% collect information about response counts and frequencies
strainCounts12 = NaN(numStrains,4);
strainFrequency12 = NaN(numStrains,4);
strainCounts24 = NaN(numStrains,4);
strainFrequency24 = NaN(numStrains,4);

% loop through each "strain" (strainListLHS combination)
for strain_id = 1:numStrains
    warn_id = all_warnall(strain_id,:) == 1;
    if sum(warn_id) > errorThresh
        CST_post = get_CST_post(squeeze(all_select_outcomes(strain_id,warn_id,time_id(1),2:end)));
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts);
        strainCounts12(strain_id,:) = Counts;
        strainFrequency12(strain_id,:) = Frequency;

        CST_post = get_CST_post(squeeze(all_select_outcomes(strain_id,warn_id,time_id(2),2:end)));
        Counts = [sum(CST_post == 1),sum(CST_post == 2),sum(CST_post == 3),sum(CST_post == 4)];
        Frequency = Counts ./ sum(Counts);
        strainCounts24(strain_id,:) = Counts;
        strainFrequency24(strain_id,:) = Frequency;
    end
end
disp('Maximum Response Frequency 12wks: nAB, Li, oLB, Probiotic:')
disp(max(strainFrequency12))

disp('Maximum Response Frequency 12wks: nAB, Li, oLB, Probiotic:')
disp(max(strainFrequency24))

%% 3. DETERMINE WHICH PROBIOTICS MATCH LACTIN-V RESULTS
ClinicalRecur = [46 59]; % Number of patients that were BV+ by week 12, week 24
ClinicalCure = [87 63]; % Number of patients that were not BV+ by week 12, week 24
ClinicalN = ClinicalCure + ClinicalRecur; % Total number of patients with known outcome
ClinicalResults = (ClinicalRecur./ClinicalN)*100; % Frequency of Recurrence (BV+/nAB dominance)

% Loop through results, collect distance from clinical result
tot_dist = NaN(size(strainFrequency12,1),1);
for i = 1:size(strainFrequency12,1)
    tmp12 = strainFrequency12(i,1)*100; % Frequency of nAB dominated states at 12 weeks
    tmp24 = strainFrequency24(i,1)*100; % Frequency of nAB dominated states at 24 weeks

    d12 = abs(tmp12-ClinicalResults(1)); % absolute distance at 12 week
    d24 = abs(tmp24-ClinicalResults(2)); % absolute distance at 24 week
    
    tot_dist(i) = d12 + d24; % total distance
end

errorThresh = 500; % threshold for the number of successful runs (244 gives 1000 "good" runs)

er_id = numODEwarnings <= errorThresh;
histogram(tot_dist(er_id))

% quick viz of results:
title('Absolute Distance from Clinical Observation')
ylabel('Frequency')

% select acceptable deviation from clinical results
dis_id = tot_dist < 3; % total percent different
LV_strains = strainListLHS(dis_id & er_id,:);

%% 4. PLOT CANDIDATE STRAIN PARAMETERS
subplot(1,4,1)
boxchart(LV_strains(:,1));
hold on
x = repmat(1,size(LV_strains,1),1);
s1 = swarmchart(x,LV_strains(:,1),[],'red','XJitterWidth',0.5);
xticklabels(param_names(PROBinfo.strainCharIndexes(1)))
ylabel('Parameter Value')
title('Growth Term')

subplot(1,4,[2,4])
boxchart(LV_strains(:,2:end))
hold on
x = repmat(1:size(LV_strains,2)-1,size(LV_strains,1),1);
s2 = swarmchart(x,LV_strains(:,2:end),[],'red','XJitterWidth',0.5);
xticklabels(param_names(PROBinfo.strainCharIndexes(2:end)))
ylabel('Parameter Value')
hold on
yline(0,'k')
title('Interaction Terms')

set(gcf,'units','pixels','position',[1 1 900 300])

%% 5. COMPARE WITH VOLCANOS
matrix1 = LV_strains;
matrix2 = strainListLHS(~dis_id & er_id,:);
alpha = 0.1;
offset = 0.04;
classes = {'LV Strains','Non-LV Strains'};
[PrismFormat,SummaryStats] = plot_Volcano(matrix1,matrix2,alpha,offset,...
    param_names(PROBinfo.strainCharIndexes),classes);
set(gca,'fontsize',12)

%% PCA
X = strainFrequency12; % or strainFrequency24
[~,lb] = max(X > 0.4,[],2);


Xblock = [];
yblock = [];
lvid = [];
c = 1;
for i = [1 2 3 4]
    tmp = strainListLHS(lb == i & numODEwarnings <= errorThresh,:);
    lvid = [lvid; tot_dist(lb == i & numODEwarnings <= errorThresh,:)];
    Xblock = [Xblock;tmp];
    yblock = [yblock;repmat(c,size(tmp,1),1)];
    c = c + 1;
end

idx = lvid < 3;

[coeff,score,latent,tsquared,explained,mu] = pca(Xblock);

pcn1 = 1;
pcn2 = 2;
xnames = param_names(PROBinfo.strainCharIndexes);
Y = yblock;
classes = SSnms;
plot_pca(score,coeff,explained,pcn1,pcn2,xnames,Y,classes)
hold on 
plot(score(idx,1),score(idx,2),'ok','MarkerSize',10)

%% 6. RUN PLS-DA
errorThresh = 500;
X = strainFrequency24; % or strainFrequency24
[~,lb] = max(X,[],2);
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
% filettl = 'XXXXXXXXXXX-12wk-LV-ResponseTypes';
filettl = 'XXXXXXXXXXX-24wk-LV-BVefficacy';
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
