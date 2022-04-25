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

% load data
load(sprintf('%sanalysis/PFanatLabels_29pat.mat', workingDir));
load(sprintf('%sanalysis/STRFmetrics_HFA_29pat.mat', workingDir));
load(sprintf('%sanalysis/STRFmat_HFA_29pat.mat', workingDir));

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


%
[~, idxBest] = sort(metrics(:, 6), 'descend');

offset = 0;
figure('Position', [193 47 1500 1500], 'Color', [1 1 1], 'DefaultAxesFontSize', 5, 'resize', 'off');
for idxTMP = 1:16
    idx = idxBest(idxTMP+offset);
    STRF = reshape(STRFmat(:, idx), 32, 75);
    subplot(4,4,idxTMP);
    imagesc(STRF, [-5 5]);
    colormap jet; axis xy; title(sprintf('%s e%i - R2=%.3g', inputArgList{metrics(idx, 1), 1}, metrics(idx, 2), metrics(idx, 6)));
end

values = nan(2668, 1);
values(idxSig(idxBest(1:20))) = 1:20;
PF_plotValOnAnat(inputArgList, 'HFA', values, 'MNI', 2, 0, 'rh');


% left candidates
candidates = [10 63; 11 67; 14 20; 26 62; 37 70];
[~, candidates(:, 1)] = ismember(candidates(:,1), cellfun(@(x) str2double(x(4:end)), inputArgList(:, 1)));
idxList = find(ismember(metrics(:, 1:2), candidates, 'rows'));

values = nan(2668, 1);
values(idxSig(idxList)) = 1:length(idxList);
PF_plotValOnAnat(inputArgList, 'HFA', values, 'MNI', 2, 0, 'lh', [nan nan], 0, 'jet');

figure('Color', [1 1 1], 'Position', [45 546 2295 362], 'DefaultAxesFontSize', 5);
for idxTMP = 1:length(idxList)
    idx = idxList(idxTMP);
    STRF = reshape(STRFmat(:, idx), 32, 75);
    subplot(1,5,idxTMP);
    imagesc(STRF, [-5 5]);
    colormap jet; axis xy; title(sprintf('%s e%i - %s - R2=%.3g', inputArgList{metrics(idx, 1), 1}, metrics(idx, 2), labels{idx}, metrics(idx, 6)));
end


% right candidates
candidates = [8 20; 9 91; 40 60; 44 15; 62 195];
[~, candidates(:, 1)] = ismember(candidates(:,1), cellfun(@(x) str2double(x(4:end)), inputArgList(:, 1)));
idxList = find(ismember(metrics(:, 1:2), candidates, 'rows'));

values = nan(2668, 1);
values(idxSig(idxList)) = 1:length(idxList);
PF_plotValOnAnat(inputArgList, 'HFA', values, 'MNI', 2, 0, 'rh', [nan nan], 0, 'jet');

figure('Color', [1 1 1], 'Position', [45 546 2295 362], 'DefaultAxesFontSize', 5);
for idxTMP = 1:length(idxList)
    idx = idxList(idxTMP);
    STRF = reshape(STRFmat(:, idx), 32, 75);
    subplot(1,5,idxTMP);
    imagesc(STRF, [-5 5]);
    colormap jet; axis xy; title(sprintf('%s e%i - %s - R2=%.3g', inputArgList{metrics(idx, 1), 1}, metrics(idx, 2), labels{idx}, metrics(idx, 6)));
end
