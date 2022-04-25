function [reconMat, r2actual, wavRecon] = PF_plotRecon(filenameRoot, paramsRecon)

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

if nargin < 2
    paramsRecon = [-1 1 0]; % threshPrctile, flagFig, forceRecon
end
paramsRecon = num2cell(paramsRecon);
[threshPrctile, flagFig, forceRecon] = deal(paramsRecon{:});

figPos = [6 479 2550 817]; % [6 701 2550 543]; %[6 569 3430 687];
denoisingExp = 2;
if forceRecon > 0
    nIterRecon = forceRecon;
else
    nIterRecon = 10;
end

load(sprintf('%s_stimuli/thewall1_stim128.mat', workingDir), 'stim128', 'CF128');
nFiles = length(CF128);
fTicks = [1 12:24:128 128];
fTickLabels = [180 250*2.^(0:4) 7200];

fileList = dir(sprintf('%srecon/%s*', workingDir, filenameRoot));
fileList = {fileList(:).name};
[actualIdx, idxSort] = sort(cellfun(@(x) str2double(x(length(filenameRoot)+1:end-4)), fileList));
fileList = fileList(idxSort);
nActualFiles = length(fileList);
if nActualFiles ~= nFiles
    fprintf('/!\\ only %i out of %i models have been computed so far /!\\\n', nActualFiles, nFiles);
end

load(sprintf('%srecon/%s', workingDir, fileList{1}), 'save_r', 'save_y_pred', 'params', 'hyperparam_grid');
nRs = length(save_r);
L = size(save_y_pred, 2);

Lstim = size(stim128, 1);
tb = linspace(0, (Lstim-1)/100, Lstim);
[~, idxTestSet] = min(abs(tb' - double(params.fixed_test)));
idxTestSet = idxTestSet(1):idxTestSet(2)-1;
if length(idxTestSet) ~= L
    idxTestSet = idxTestSet(1:L);
end
stimTarget = stim128(idxTestSet, :)';
tbTestSet = tb(idxTestSet);

hpFields = {'learning_rate', 'l1_ratio', 'hidden_layer_sizes', 'alpha'};
hpScales = {'log', 'linear', 'linear', 'log'};
if exist('hyperparam_grid', 'var')
    hpField = fieldnames(hyperparam_grid);
    hpScale = hpScales(ismember(hpFields, hpField));
else
    switch params.algo
        case 'rMLRwES'
            hpField = {'learning_rate'};
            hpScale = {'log'};
        case {'ridge', 'lasso'}
            hpField = {'alpha'};
            hpScale = {'linear'};
        case 'elasticNet'
            hpField = {'l1_ratio', 'alpha'};
            hpScale = {'linear'};
        case 'MLP'
            hpField = {'hidden_layer_sizes', 'alpha'};
            hpScale = {'linear', 'log'};
    end
end
nHp = length(hpField);

r2Mat = nan(nFiles, nRs);
hpCell = cell(nFiles, nHp);
yPredMat = nan(nFiles, nRs, L);
for idxFile = 1:nActualFiles
    load(sprintf('%srecon/%s', workingDir, fileList{idxFile}), 'save_r2', 'params', 'save_y_pred');
    r2Mat(actualIdx(idxFile), :) = save_r2;
    for idxHp = 1:nHp
        hpCell{actualIdx(idxFile), idxHp} = params.(hpField{idxHp});
    end
    yPredMat(actualIdx(idxFile), :, :) = save_y_pred;
end

reconMat = zeros(nFiles, L);
save_nRs = nan(nFiles, 1);
for idxFile = 1:nActualFiles
    idxTMP = actualIdx(idxFile);
    metricsTMP = r2Mat(idxTMP, :);
    if threshPrctile > 0
        idxBest = metricsTMP > prctile(metricsTMP, threshPrctile);
        threshStr = sprintf(' - best %.03g out of %.03g models', sum(idxBest), nRs);
        save_nRs(idxTMP) = sum(idxBest);
    elseif threshPrctile == 0
        idxBest = true(1, nRs);
        threshStr = sprintf(' - all %.03g models', nRs);
        save_nRs(idxTMP) = sum(idxBest);
    else
        [~, idxSort] = sort(metricsTMP);
        effectiveMetricsTMP = zeros(nRs, 1);
        for idxRs = 1:nRs
            effectiveMetricsTMP(idxRs) = computeR2(stimTarget(idxTMP, :), squeeze(mean(yPredMat(idxTMP, idxSort(idxRs:end), :), 2)));
        end
        [~, idxMax] = max(effectiveMetricsTMP);
        idxBest = false(1, nRs);
        idxBest(idxSort(idxMax:end)) = true;
        save_nRs(idxTMP) = nRs - idxMax + 1;
        threshStr = ' - adaptive nRs';
    end
    reconMat(idxTMP, :) = squeeze(mean(yPredMat(idxTMP, idxBest, :), 2));
end

r2actual = nan(nFiles, 1);
for idxF = 1:nActualFiles
    r2actual(actualIdx(idxF)) = computeR2(stimTarget(actualIdx(idxF), :), reconMat(actualIdx(idxF), :));
end
r2d = corr2(stimTarget(actualIdx, :), reconMat(actualIdx, :));

idxEmpty = find(isnan(r2actual));
if ~isempty(idxEmpty)
    idxBreak = [0 find(diff(idxEmpty) > 1)' length(idxEmpty)];
    for idxB = 1:length(idxBreak)-1
        idxEmptyTMP = idxEmpty(idxBreak(idxB)+1:idxBreak(idxB+1));
        reconMat(idxEmptyTMP, :) = .8.*repmat(1.*mean(reconMat([idxEmptyTMP(1)-1 min([idxEmptyTMP(end)+1 nFiles])], :)), length(idxEmptyTMP), 1);
    end
end


if flagFig > 0
    mT = [mean(stimTarget(:)) median(stimTarget(:)) std(stimTarget(:)) min(stimTarget(:)) max(stimTarget(:))];
    mR = [mean(reconMat(:), 'omitnan') median(reconMat(:), 'omitnan') std(reconMat(:), 'omitnan') min(reconMat(:)) max(reconMat(:))];
    
    if ~exist('hyperparam_grid', 'var')
        if nHp == 2
            hpStr = sprintf(' - %s: %s - %s: %g', hpField{1}, num2str(params.(hpField{1})), hpField{2}, params.(hpField{2}));
        else
            hpStr = sprintf(' - %s: %g', hpField{1}, params.(hpField{1}));
        end
        nHp = 0;
    else
        hpStr = '';
    end
    hFig = figure('Position', figPos, 'Color', [1 1 1], 'DefaultAxesFontSize', 5, 'Name', filenameRoot);
    
    subplot(2,5+nHp,1:3); imagesc(tbTestSet, 1:nFiles, stimTarget, [mT(4) mT(1)+2*mT(3)]); axis xy; ylim([1 nFiles]); set(gca, 'YTick', fTicks, 'YTickLabel', fTickLabels); ylabel('frequency (Hz)'); colorbar; title('original spectrogram'); xlabel('time (s)');
    titleStr = sprintf('reconstructed spectrogram%s - corr2=%.3g%s', hpStr, r2d, threshStr);
    subplot(2,5+nHp,(1:3)+5+nHp); imagesc(tbTestSet, 1:nFiles, reconMat, [mR(4) mR(1)+2*mR(3)]); axis xy; colorbar; ylim([1 nFiles]); set(gca, 'YTick', fTicks, 'YTickLabel', fTickLabels); title(titleStr, 'Interpreter', 'none');
    uicontrol('parent', hFig, 'style', 'push', 'units', 'pix', 'position', [10 250 120 40], 'fontsize', 5, 'string', 'wav orig', 'callback', @pb1t_callback);
    uicontrol('parent', hFig, 'style', 'push', 'units', 'pix', 'position', [10 130 120 40], 'fontsize', 5, 'string', 'wav recon', 'callback', @pb1r_callback);
    uicontrol('parent', hFig, 'style', 'push', 'units', 'pix', 'position', [10 210 120 40], 'fontsize', 5, 'string', 'play orig', 'callback', @pb2t_callback);
    uicontrol('parent', hFig, 'style', 'push', 'units', 'pix', 'position', [10 90 120 40], 'fontsize', 5, 'string', 'play recon', 'callback', @pb2r_callback);
    uicontrol('parent', hFig, 'style', 'push', 'units', 'pix', 'position', [10 10 120 40], 'fontsize', 5, 'string', 'stop', 'callback', @pb3_callback);
    uicontrol('parent', hFig, 'style', 'edit', 'units', 'pix', 'position', [140 210 120 40], 'fontsize', 5, 'callback', @pb5_callback);
    hText = uicontrol('parent', hFig, 'style', 'text', 'units', 'pix', 'position', [140 250 120 40], 'fontsize', 5, 'string', sprintf('nIter = %i change below', nIterRecon));
    uicontrol('parent', hFig, 'style', 'edit', 'units', 'pix', 'position', [140 90 120 40], 'fontsize', 5, 'callback', @pb4_callback);
    hText2 = uicontrol('parent', hFig, 'style', 'text', 'units', 'pix', 'position', [140 130 120 40], 'fontsize', 5, 'string', sprintf('DE = %.3g change below', denoisingExp));
   
    h = zeros(1+nHp);
    h(1) = subplot(2,5+nHp,[4 9+nHp]); plot(0:nFiles+1, zeros(nFiles + 2, 1), ':k'); hold on; plot(r2actual, '.r'); hold off; ylabel('effective r-squared'); box off; view([90 -90]); set(gca, 'XTick', fTicks, 'XTickLabel', fTickLabels); ylim([-1 1]);
    title(sprintf('mean = %.3g - median = %.3g\nmin = %.3g - max = %.3g', mean(r2actual, 'omitnan'), median(r2actual, 'omitnan'), min(r2actual), max(r2actual)));
    h(2) = subplot(2,5+nHp,[5 10+nHp]); plot(save_nRs, '.b'); ylabel('nRs'); box off; view([90 -90]); set(gca, 'XTick', fTicks, 'XTickLabel', []); ylim([0 nRs+1]);
    title(sprintf('mean = %.3g - median = %.3g\nmin = %.3g - max = %.3g', mean(save_nRs, 'omitnan'), median(save_nRs, 'omitnan'), min(save_nRs), max(save_nRs)));
    
    if nHp == 1
        hpMat = nan(nFiles, 1);
        idxOk = find(~cellfun(@isempty, hpCell));
        hpMat(idxOk) = cell2mat(hpCell(idxOk));
        h(3) = subplot(2,5+nHp, [6 11+nHp]); plot(hpMat, '.k'); ylabel(hpField{1}, 'Interpreter', 'none'); set(gca, 'YScale', hpScale{1});
        box off; view([90 -90]); set(gca, 'XTick', fTicks, 'XTickLabel', []); yticks(hyperparam_grid.(hpField{1})); title('hyperparameter tuning');
    elseif nHp > 1
        for idxHp = 1:nHp
            h(2+idxHp) = subplot(2,5+nHp,[5+idxHp 10+nHp+idxHp]);
            hpTMP = hpCell(:, idxHp);
            try
                uniqueLabel = cellfun(@num2str, hyperparam_grid.(hpField{idxHp}), 'un', 0);
            catch
                uniqueLabel = sprintfc('%.03g', hyperparam_grid.(hpField{idxHp}));
            end
            hpTMPstr = cellfun(@num2str, hpTMP, 'un', 0);
            hpTMPidx = cellfun(@(x) find(ismember(uniqueLabel, x)), hpTMPstr, 'un', 0);
            hpTMPidx2 = nan(nFiles, 1);
            idxOk = find(~cellfun(@isempty, hpTMPidx));
            hpTMPidx2(idxOk) = cell2mat(hpTMPidx(idxOk));
            plot(hpTMPidx2, '.k'); hold on;
            set(gca, 'YScale', hpScale{idxHp});
            ylabel(hpField{idxHp}, 'Interpreter', 'none');
            box off; view([90 -90]); title('hyperparameter tuning'); set(gca, 'XTick', fTicks, 'XTickLabel', []);
            ylim([0 length(uniqueLabel)+1]); set(gca, 'YTick', 1:length(uniqueLabel), 'YTickLabel', uniqueLabel);
        end
    end
    linkaxes(h, 'x'); xlim([1 nFiles]);
    
    if forceRecon > 0
        [~, wavRecon] = PF_aud2wav(reconMat'.^denoisingExp, 0, forceRecon);
    else
        wavRecon = [];
    end
    
    S = struct('stimTarget', stimTarget, 'stimRecon', reconMat, 'nIterRecon', nIterRecon, 'denoisingExp', denoisingExp);
    if exist('wavRecon', 'var')
        S.wavRecon = wavRecon;
    end
    S.hText = hText;
    S.hText2 = hText2;
    guidata(hFig, S);
end

if forceRecon > 0 && ~exist('wavRecon', 'var')
    [~, wavRecon] = PF_aud2wav(reconMat'.^denoisingExp, 0, forceRecon);
else
    wavRecon = [];
end


function pb1t_callback(hObject, eventdata)
fprintf('Computing wav from target\n');
data = guidata(hObject);
[~, wavTarget] = PF_aud2wav(data.stimTarget', 0, data.nIterRecon);
data.wavTarget = wavTarget;
guidata(hObject, data);
pb2t_callback(hObject);


function pb1r_callback(hObject, eventdata)
fprintf('Computing wav from recon\n');
data = guidata(hObject);
[~, wavRecon] = PF_aud2wav(data.stimRecon'.^data.denoisingExp, 0, data.nIterRecon);
data.wavRecon = wavRecon;
guidata(hObject, data);
pb2r_callback(hObject);


function pb2t_callback(hObject, eventdata)
fprintf('Playing target wav\n');
data = guidata(hObject);
soundsc(data.wavTarget, 16000)


function pb2r_callback(hObject, eventdata)
fprintf('Playing recon wav\n');
data = guidata(hObject);
soundsc(data.wavRecon, 16000)


function pb3_callback(hObject, eventdata)
clear sound

function pb4_callback(hObject, eventdata)
fprintf('Updating denoisingExponent value\n');
data = guidata(hObject);
data.denoisingExp = str2double(hObject.String);
data.hText2.String = sprintf('DE = %s change below', hObject.String);
hObject.String = '';
guidata(hObject, data);

function pb5_callback(hObject, eventdata)
fprintf('Updating denoisingExponent value\n');
data = guidata(hObject);
data.nIterRecon = str2double(hObject.String);
data.hText.String = sprintf('nIter = %i change below', data.nIterRecon);
hObject.String = '';
guidata(hObject, data);
