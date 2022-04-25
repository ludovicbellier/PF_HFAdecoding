function [R2out, rOut] = PF_plotDecodingMetrics(patientCode, suffix, modelType)

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

flagFig = 0;

metrics = PF_getModelMetrics(patientCode, suffix, modelType);
nFBins = size(metrics, 1);
idxNonSigR = find(metrics(:, 3) == 0);
idxNonSigR2 = find(metrics(:, 4) == 0);

lim_nan = 200;
if any(metrics(:, 7) > lim_nan) % display message if model with nan_lim exceeded
    fprintf('/!\\ %s_%s has %i models with more than %i failed iterations\n', patientCode, suffix, sum(metrics(:, 7) > 200), lim_nan);
end

idx = 1;
while ~exist('nRs', 'var')
    try
        load(sprintf('%s%s/%s_%s_f%i.mat', workingDir, modelType, patientCode, suffix, idx), 'save_r');
        nRs = length(save_r);
    catch
        idx = idx + 1;
    end
end
superSaveR = zeros(nRs, nFBins);
superSaveR2 = zeros(nRs, nFBins);
for idx = 1:nFBins
    try
        load(sprintf('%s%s/%s_%s_f%i.mat', workingDir, modelType, patientCode, suffix, idx), 'save_r', 'save_r2');
    catch
        save_r = nan(nRs, 1);
        save_r2 = nan(nRs, 1);
    end
    superSaveR(:, idx) = save_r;
    superSaveR2(:, idx) = save_r2;
end

if nFBins == 128
    idxXTick = [1 16:16:128];
else
    idxXTick = [1 4:4:32];
end
if flagFig > 0
    figure('Position', [6 -204 1200 884], 'Name', sprintf('%s/%s_%s', modelType, patientCode, suffix), 'DefaultAxesFontSize', 5);
    subplot(211); hold on;
    for idx = -.8:.2:.8
        line([0 nFBins+1], [idx idx], 'Color', [.9 .9 .9]);
    end
    line([0 nFBins+1], [0 0], 'Color', [.7 .7 .7]);
    plot(metrics(:, 1), 'r', 'linewidth', 2);
    plot(idxNonSigR, metrics(idxNonSigR, 1), '+b', 'markersize', 10);
    boxplot(superSaveR);
    xlim([0 nFBins+1]); ylim([-1 1]);
    set(gca, 'XTick', idxXTick, 'XTickLabel', idxXTick, 'box', 'off');
    title(sprintf('mean r=%.2g - %i sig fBins', mean(metrics(:, 1)), sum(metrics(:, 3))));
    
    subplot(212); hold on;
    for idx = -.8:.2:.8
        line([0 nFBins+1], [idx idx], 'Color', [.9 .9 .9]);
    end
    line([0 nFBins+1], [0 0], 'Color', [.7 .7 .7]);
    plot(metrics(:, 2), 'r', 'linewidth', 2);
    plot(idxNonSigR2, metrics(idxNonSigR2, 2), '+b', 'markersize', 10);
    boxplot(superSaveR2);
    xlim([0 nFBins+1]); ylim([-1 1]);
    set(gca, 'XTick', idxXTick, 'XTickLabel', idxXTick, 'box', 'off');
    title(sprintf('mean R2=%.2g - %i sig fBins', mean(metrics(:, 2)), sum(metrics(:, 4))));
end
rOut = metrics(:, 1);
R2out = metrics(:, 2);