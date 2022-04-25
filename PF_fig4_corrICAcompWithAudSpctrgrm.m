global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

load(sprintf('%sanalysis/ICAcomps.mat', workingDir), 'comp_ICA_2D');
[nFreq, nLags, nComps] = size(comp_ICA_2D);

load(sprintf('%s_stimuli/thewall1_stim32.mat', workingDir));
L2 = size(stim32, 1) - nLags + 1;

compCorr = zeros(L2, nComps);
for idxComp = 1:nComps
    for idx = 1:L2
        STMP = stim32(idx+(0:nLags-1), :);
        CTMP = squeeze(comp_ICA_2D(:, :, idxComp))';
        compCorr(idx, idxComp) = corr2(STMP, CTMP);
    end
end
compCorr = [nan(nLags-1, nComps); compCorr];

%% Fig. 4A and C
thresh = 0;
offset = 0;
compNames = {'onset', 'sustained', 'late onset'};
load(sprintf('%s_stimuli/thewall1_stim128.mat', workingDir));
L = size(stim128, 1);
tb = linspace(0, (L-1)/100, L);
h = zeros(nComps+1, 1);
figure('Position', [72 140 1837 1116], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
h(1) = subplot(nComps+2,1,1:2); imagesc(tb, 1:128, log10(stim128)', [-2 0]); axis xy;
yTick = [1 12 36 60 84 108 128];
set(gca, 'YTick', yTick, 'YTickLabel', round(CF128(yTick)));
saveYLim = zeros(3, 2);
for idxComp = 1:nComps
    h(1+idxComp) = subplot(nComps+2,1,2+idxComp);
    plot(tb, compCorr(:, idxComp+offset)); hold on;
    plot(tb(compCorr(:, idxComp+offset)>thresh), compCorr(compCorr(:, idxComp+offset)>thresh, idxComp+offset), 'r.', 'MarkerSize', 2);
    box off; title(compNames{idxComp});
    saveYLim(idxComp, :) = get(h(1+idxComp), 'YLim');
end
linkaxes(h, 'x'); xlim(tb([1 end]));

% for idx = 1:19
%     xlim([0 10]+10*(idx-1));
%     input('');
% end

%% Fig. 4D
set(gcf, 'Position', [72 140 573 1116]);
xlim([32 42]);
for idxComp = 1:3
    set(h(1+idxComp), 'YLim', saveYLim(idxComp, :));
end

%% Fig. 4E
xlim([100 110]);
for idxComp = 1:3
    set(h(1+idxComp), 'YLim', saveYLim(idxComp, :));
end

%% Fig. 4B
figure('Position', [72 338 274 918], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
for i = 1:nComps
    subplot(nComps,1,i);
    imagesc(squeeze(comp_ICA_2D(:,:,i)), [-25 25]); axis xy;
    set(gca, 'XTick', [], 'YTick', []);
    colormap jet;
end
