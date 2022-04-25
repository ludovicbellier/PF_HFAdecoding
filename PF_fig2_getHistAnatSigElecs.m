global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

inputArgList = {'AMC006', 1; 'AMC007', 1; 'AMC008', 1; 'AMC009', 1; 'AMC010', 1;...
'AMC011', 1; 'AMC012', 2; 'AMC013', 1; 'AMC014', 1; 'AMC015', 1;...
'AMC017', 1; 'AMC018', 2; 'AMC019', 1; 'AMC022', 1; 'AMC026', 1;...
'AMC027', 1; 'AMC028', 1; 'AMC029', 3; 'AMC031', 1; 'AMC032', 1;...
'AMC033', 1; 'AMC037', 1; 'AMC038', 1; 'AMC039', 2; 'AMC040', 2;...
'AMC041', 2; 'AMC044', 1; 'AMC045', 1; 'AMC062', 1};

load(sprintf('%sanalysis/PFanatLabels_29pat.mat', workingDir));
load(sprintf('%sanalysis/STRFmetrics_HFA_29pat.mat', workingDir));


%% 2A - Coverage
PF_plotValOnAnat(inputArgList, 'HFA', '', 'mni', 0, 0, 'lh');
PF_plotValOnAnat(inputArgList, 'HFA', '', 'mni', 0, 0, 'rh');


%% Set anatomical regions
idxColSigR = ismember(metricsCols, 'sigR');
idxColR2 = ismember(metricsCols, 'R2');
threshR2 = 0;
if ~isnan(threshR2)
    idxSig = find(metrics(:, idxColSigR) == 1 & metrics(:, idxColR2) > threshR2 & sum(metrics(:, end-1:end), 2) == 0);
else
    idxSig = find(metrics(:, idxColSigR) == 1 & sum(metrics(:, end-1:end), 2) == 0);
end
labels = labels(idxSig);

metrics = metrics(idxSig, [1:9 16:18]);
metricsCols = metricsCols([1:9 16:18]);

idx_r_STG = find(ismember(labels, 'ctx-rh-superiortemporal'));
idx_r_sensorimotor = find(ismember(labels, {'ctx-rh-precentral', 'ctx-rh-postcentral'}));
idx_r_IFG = find(ismember(labels, {'ctx-rh-parsopercularis', 'ctx-rh-parstriangularis', 'ctx-rh-parsorbitalis'}));
idx_r_supramarginal = find(ismember(labels, 'ctx-rh-supramarginal'));
idx_r_otherTemporal = find(ismember(labels, 'ctx-rh-middletemporal'));
idx_r_otherFrontal = find(ismember(labels, {'ctx-rh-rostralmiddlefrontal', 'ctx-rh-superiorfrontal', 'ctx-rh-caudalmiddlefrontal'}));
idx_r_list = {idx_r_STG, idx_r_sensorimotor, idx_r_IFG, idx_r_supramarginal, idx_r_otherTemporal, idx_r_otherFrontal};
C_r = cellfun(@length, idx_r_list);

idx_l_STG = find(ismember(labels, 'ctx-lh-superiortemporal'));
idx_l_sensorimotor = find(ismember(labels, {'ctx-lh-precentral', 'ctx-lh-postcentral'}));
idx_l_IFG = find(ismember(labels, {'ctx-lh-parsopercularis', 'ctx-lh-parstriangularis', 'ctx-lh-parsorbitalis'}));
idx_l_supramarginal = find(ismember(labels, 'ctx-lh-supramarginal'));
idx_l_otherTemporal = find(ismember(labels, {'ctx-lh-bankssts', 'ctx-lh-inferiortemporal', 'ctx-lh-middletemporal'}));
idx_l_otherFrontal = find(ismember(labels, {'ctx-lh-rostralmiddlefrontal', 'ctx-lh-superiorfrontal', 'ctx-lh-caudalmiddlefrontal'}));
idx_l_list = {idx_l_STG, idx_l_sensorimotor, idx_l_IFG, idx_l_supramarginal, idx_l_otherTemporal, idx_l_otherFrontal};
C_l = cellfun(@length, idx_l_list);


%% 2B - Plot significant electrodes with marker size mapped to r
allIdx = {idx_l_STG; idx_l_sensorimotor; idx_l_IFG; idx_l_supramarginal; idx_l_otherTemporal; idx_l_otherFrontal;
    idx_r_STG; idx_r_sensorimotor; idx_r_IFG; idx_r_supramarginal; idx_r_otherTemporal; idx_r_otherFrontal};
Ntotal = 2668;
anatValues = zeros(Ntotal, 1);
for idx = 1:length(allIdx)
    anatValues(idxSig(allIdx{idx})) = idx;
end

colorList = [213 32 30; 1 147 74; 2 117 178; 219 131 26; 235 210 1; 105 50 127]./255;
colorList = colorList./max(colorList, [], 'all');

sizeRange = [50 300];

[~, idxOrig] = PF_plotValOnAnat(inputArgList, 'HFA', anatValues, 'mni', 2, 0, 'lh', [1 6], 0, colorList);
hScatter = findobj('type', 'Scatter');
idxLeft = ismember(idxSig, idxOrig);
sizeValues = metrics(idxLeft, 5);
set(hScatter, 'SizeData', rescale(sizeValues, sizeRange(1), sizeRange(2)));

[~, idxOrig] = PF_plotValOnAnat(inputArgList, 'HFA', anatValues, 'mni', 2, 0, 'rh', [7 12], 0, colorList);
hScatter = findobj('type', 'Scatter');
idxRight = ismember(idxSig, idxOrig);
sizeValues = metrics(idxRight, 5);
set(hScatter, 'SizeData', rescale(sizeValues, sizeRange(1), sizeRange(2)));


%% 2C - pie chart with separate l and r
figure('Position', [67 638 2439 582], 'Color', [1 1 1], 'DefaultAxesFontSize', 4);
newC = [C_l; C_r];
newC = newC(:);
nC = length(newC);
Clabels = cell(nC, 1);
for i = 1:nC
    Clabels{i} = sprintf('%i', newC(i));
end
p = pie(newC, zeros(nC, 1), Clabels);
colorList2 = repelem(colorList, 2, 1);
colorList2(2:2:end) = colorList2(2:2:end).*.7;
for idxC = 1:nC
    p(2*idxC-1).FaceColor = colorList2(idxC, :);
end
vals = get(gca, 'Children');
for i = 1:2:24
    set(vals(i), 'FontSize', 20);
end


%% Pool together "other" elecs
idx_r_all = find(metrics(:, 3) == 2);
idx_r_other = setdiff(idx_r_all, [idx_r_STG; idx_r_sensorimotor; idx_r_IFG]);
idx_r_list = {idx_r_STG, idx_r_sensorimotor, idx_r_IFG, idx_r_other};
C_r = cellfun(@length, idx_r_list);
r_mean_r = cellfun(@(x) mean(metrics(x, 5)), idx_r_list);
r_err_r = cellfun(@(x) std(metrics(x, 5))/sqrt(length(x)), idx_r_list);
R2_mean_r = cellfun(@(x) mean(metrics(x, 6)), idx_r_list);
R2_err_r = cellfun(@(x) std(metrics(x, 6))/sqrt(length(x)), idx_r_list);

idx_l_all = find(metrics(:, 3) == 1);
idx_l_other = setdiff(idx_l_all, [idx_l_STG; idx_l_sensorimotor; idx_l_IFG]);
idx_l_list = {idx_l_STG, idx_l_sensorimotor, idx_l_IFG, idx_l_other};
C_l = cellfun(@length, idx_l_list);
r_mean_l = cellfun(@(x) mean(metrics(x, 5)), idx_l_list);
r_err_l = cellfun(@(x) std(metrics(x, 5))/sqrt(length(x)), idx_l_list);
R2_mean_l = cellfun(@(x) mean(metrics(x, 6)), idx_l_list);
R2_err_l = cellfun(@(x) std(metrics(x, 6))/sqrt(length(x)), idx_l_list);


%% 2D - bar plot for r values
colorList3 = colorList2;
colorList3(7:end, :) = [];
colorList3(7:8, :) = [1 1 1; .7 .7 .7];
figure('Position', [625 490 1460 740], 'Color', [1 1 1], 'DefaultAxesFontSize', 4);
nC = length(C_r);
hB = zeros(2*nC, 1);
rVals = [r_mean_l; r_mean_r];
rVals = rVals(:);
offsetVal = .21;
offsets = [-1; 1]*offsetVal+(1:nC);
offsets = offsets(:);
for i = 1:nC*2
    hB(i) = bar(offsets(i), rVals(i), 'BarWidth', .4, 'FaceColor', colorList3(i, :)); hold on;
end
hold on; errorbar((1:nC)-offsetVal, r_mean_l, [], r_err_l, 'LineStyle', 'none', 'Color', [0 0 0]);%colorList(1, :));
errorbar((1:nC)+offsetVal, r_mean_r, [], r_err_r, 'LineStyle', 'none', 'Color', [0 0 0]);% colors(2, :));
set(gca, 'XTick', 1:4, 'XTickLabels', {'STG', 'SMC', 'IFG', 'other'});
ylabel('Pearson''s r'); box off; xlim([.5 nC+.5]);
set(gca, 'YGrid', 'on');


%% Stats - anovan / unbalanced two-way anova with all labels
idxLeft = vertcat(idx_l_list{:});
idxRight = vertcat(idx_r_list{:});
y = [atanh(metrics(idxLeft, 5)); atanh(metrics(idxRight, 5))];
g1 = [ones(length(idxLeft), 1); 2*ones(length(idxRight), 1)];
g2LeftIdx = [0 cumsum(C_l)];
g2RightIdx = [0 cumsum(C_r)];
g2Left = zeros(g2LeftIdx(end), 1);
g2Right = zeros(g2RightIdx(end), 1);
for i = 1:length(C_l)
    g2Left(g2LeftIdx(i)+1:g2LeftIdx(i+1)) = i;
    g2Right(g2RightIdx(i)+1:g2RightIdx(i+1)) = i;
end
g2 = [g2Left; g2Right];
g1Names = {'Left', 'Right'};
g1 = g1Names(g1);
g2Names = {'STG', 'SMC', 'IFG', 'other'};%'supramarginal', 'otherTemporal', 'otherFrontal'};
g2 = g2Names(g2);
[p, tbl, stats] = anovan(y, {g1 g2}, 'model', 'linear', 'varnames', {'laterality', 'area'});
% csvwrite('PF_fig1_rValuesGivenAreasAndLat.csv', [g1 g2 y]);
[c, m] = multcompare(stats);
[c, m] = multcompare(stats, 'dimension', 2);
