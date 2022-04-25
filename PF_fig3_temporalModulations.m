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

load(sprintf('%sanalysis/STRFmetrics_HFA_29pat.mat', workingDir));
idxColSigR = ismember(metricsCols, 'sigR');
idxColR2 = ismember(metricsCols, 'R2');
threshR2 = 0;
if ~isnan(threshR2)
    idxSig = find(metrics(:, idxColSigR) == 1 & metrics(:, idxColR2) > threshR2 & sum(metrics(:, end-1:end), 2) == 0);
else
    idxSig = find(metrics(:, idxColSigR) == 1 & sum(metrics(:, end-1:end), 2) == 0);
end
metrics = metrics(idxSig, :);

addpath('/home/knight/lbellier/DataWorkspace/_tools/nsltools/');
N = length(idxSig);
flagFig = 0;
rv = 2.^(0:.2:5);
R = length(rv);
thresh = .3;
saveIndices = zeros(N, 4);
rfSave = cell(N, 1);
idxFOI = 1:32;
idxROI = find(rv<1000/150, 1, 'last'):find(rv<1000/150, 1, 'last')+2;
load(sprintf('%sanalysis/STRFmat_HFA_29pat.mat', workingDir));
for idx = 1:N
    STRFtmp = reshape(STRFmat(:, idx), 32, 75);
    rtf = aud2tf(STRFtmp', rv, [], 100, 24, 1);
    rf = squeeze(mean(abs(rtf(R+1:end,:,:)), 2))';
    rfSave{idx} = rf;
    indexTMP = max(rf(idxFOI, idxROI), [], 'all');
    indexTMP2 = mean(rf(idxFOI, idxROI), 'all');
    indexTMP3 = max(mean(rf(idxFOI, find(rv>4, 1, 'first'):end)));
    indexTMP4 = mean(rf(1:15, idxROI), 'all');
    if flagFig > 0
        figure(2);
        subplot(221); imagesc(STRFtmp, [-5 5]); axis xy; colormap jet;
        tmp = metrics(idx, 1:2);
        title(sprintf('%s e%i', inputArgList{tmp(1), 1}, tmp(2)));
        subplot(222); imagesc(rf); axis xy
        set(gca, 'XTickLabel', round(rv(get(gca, 'XTick'))));
        title(sprintf('%.3g', indexTMP));
        subplot(223); bar(indexTMP3); ylim([0 1]);
        if indexTMP >= thresh
            title('rhythmic', 'color', 'r');
        end
        subplot(224); plot(mean(rf)); %plot(rf');%
        xlim([1 R]);
        set(gca, 'XTickLabel', round(rv(get(gca, 'XTick'))));
        title(sprintf('%.3g', indexTMP3));
        input('');
    end
    saveIndices(idx, :) = [indexTMP indexTMP2 indexTMP3 indexTMP4];
end

% save functional sets of electrodes for the ablation analysis
load(sprintf('%sanalysis/STRFcoeffICA.mat', workingDir));
idxFuncSets = coeff_ICA > 0;
idxFuncSets(:, 4) = saveIndices(:, 1) >= thresh;
funcSetNames = {'onset', 'sustained', 'lateOnset', 'rhythmic'};
save(sprintf('%sanalysis/STRFfuncSets.mat', workingDir), 'idxFuncSets', 'funcSetNames');

% Fig. 3D - example RF
figure('Position', [72 538 880 718], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
imagesc(rfSave{331}); axis xy; colormap jet;
set(gca, 'XTick', 1:5:26, 'XTickLabel', rv(1:5:26));
xlabel('Rate (Hz)'); ylabel('Frequency (Hz)'); colorbar;

% Fig. 3E - rhythmic indices plotted on the MNI brain
cMapRes = 512;
cMap = redblue(cMapRes);
cMap = [[0 0 0]; cMap(round(cMapRes/2)+1:end, :)];
values = nan(2668, 1);
values(idxSig) = (saveIndices(:, 1) - thresh).^.5;
values(values < 0) = 0;
PF_plotValOnAnat(inputArgList, 'HFA', values, 'MNI', 1, 0, 'lh', [0 nan], 0, cMap);
PF_plotValOnAnat(inputArgList, 'HFA', values, 'MNI', 1, 0, 'rh', [0 nan], 0, cMap);
