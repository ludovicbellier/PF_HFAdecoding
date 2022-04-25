global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

load(sprintf('%s_stimuli/thewall1_stim32.mat', workingDir));

suffList = {'', 'NoL', 'NoR',...
            'NoSTG', 'NoLSTG', 'NoRSTG',...
            'NoSMC', 'NoLSMC', 'NoRSMC',...
            'NoIFG', 'NoLIFG', 'NoRIFG', 'NoOther',...
            'NoOnset', 'NoLOnset', 'NoROnset',...
            'NoSustained', 'NoLSustained', 'NoRSustained',...
            'NoLateOnset', 'NoLLateOnset', 'NoRLateOnset',...
            'NoRhythmic', 'NoLRhythmic', 'NoRRhythmic'};
nSuf = length(suffList);
idxAnat = 2:13;
idxFunc = 14:25;
idxCell = {idxAnat, idxFunc};

fname = sprintf('%sanalysis/ABLATIONmetrics_HFA_%igroups.mat', workingDir, nSuf);
if exist(fname, 'file') == 0
    R2outSave = zeros(32, nSuf);
    rOutSave = zeros(32, nSuf);
    NelecsSave = zeros(nSuf, 1);
    for i = 1:nSuf
        [R2outSave(:, i), rOutSave(:, i)] = PF_plotDecodingMetrics(sprintf('supergrid%s', suffList{i}), 'TheWall1_run1_HFA', 'decoding');
        dataTMP = load(sprintf('_preprocessed_ECoG/supergrid%s_TheWall1_run1_preprocessed_HFA.mat', suffList{i}), 'ecog');
        NelecsSave(i) = size(dataTMP.ecog, 2);
    end
    save(fname, 'R2outSave', 'rOutSave', 'NelecsSave');
else
    load(fname, 'R2outSave', 'rOutSave', 'NelecsSave');
end

N = 347;
values = rOutSave - rOutSave(:, 1);
suffListTMP = cellfun(@(x, y) [x ' (-' y ')'], suffList(2:end), sprintfc('%i', N-NelecsSave(2:end))', 'un', 0);


%% Statistics, using a one-way anova
suffList2 = suffList;
suffList2{1} = 'orig';
y = atanh(values(:));
g1 = repelem(1:nSuf, 1, 32);
g1 = suffList2(g1);
[p, tbl, stats] = anovan(y, {g1}, 'model', 'linear');
c = multcompare(stats);
sigIdx = zeros(24, 1);
sigIdx(c(c(:, 1) == 1, 6) < .05) = 1;


%% Fig. 5A and B
colorList = [213 32 30; 1 147 74; 2 117 178]./255;
colorList = colorList./max(colorList, [], 'all');
colorList = repelem(colorList, 2, 1);
colorList(2:2:end) = colorList(2:2:end).*.7;
colorList(end+1:end+3, :) = [.7 .7 .7; .3 .3 .3; 1 1 1];

positions = [1:2 4:6 8:10 12:14 16 18:20 22:24 26:28 30:32];

idxColorsAnat = [7 8 1 1 2 3 3 4 5 5 6 9];
idxColorsFunc = repmat([9 7 8], 1, 4);
idxColors = [idxColorsAnat idxColorsFunc];

for idxType = 1:2
    idxTMP = idxCell{idxType} - 1;
    
    idxColorsTMP = idxColors(idxTMP);
    positionsTMP = positions(idxTMP);
    positionsTMP = positionsTMP - positionsTMP(1) + 1;
    sigTMP = find(sigIdx(idxTMP));
    
    figure('Color', [1 1 1], 'Position', [350 582 820 620], 'DefaultAxesFontSize', 5);
    line([0 max(positionsTMP)+1], [0 0], 'color', [1 1 1]*.2, 'linestyle', ':'); hold on;
    boxplot(values(:, idxTMP+1), 'Positions', positionsTMP, 'notch', 'on', 'OutlierSize', 2, 'Symbol', 'o', 'colors', 'k');%, 'colors', colorList(idxColorsTMP, :));
    box off;
    ylim([-.25 .05]); xlim([0 max(positionsTMP)+1]);
    set(gca, 'XTickLabel', suffListTMP(idxTMP-idxTMP(1)+1), 'XTickLabelRotation', 45);
    ylabel('change in Pearson''s r'); xlabel('ablation');
    
    colors2 = flipud(colorList(idxColorsTMP, :));
    h2 = findobj(gca,'Tag','Box');
    for idx = 1:length(idxColorsTMP)
        patch(get(h2(idx),'XData'), get(h2(idx),'YData'), colors2(idx, :), 'LineStyle', 'none');%, 'FaceAlpha', alphaLevels(alpha2(idx)));
    end
    
    h = get(gca, 'Children');
    h = get(h(end-1), 'Children');
    h = h(1:length(idxColorsTMP));
    for idx = 1:length(idxColorsTMP)
        set(h(idx), 'MarkerEdgeColor', colors2(idx, :));
    end
    
    set(gca,'children', flipud(get(gca,'children')));
    
    for idx = 1:length(sigTMP)
        text(positionsTMP(sigTMP(idx)), .05, '*', 'Color', 'r');
    end
    
    set(gca, 'YGrid', 'on');
end
