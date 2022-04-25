function PF_buildSTRFmat(inputArgList, suffix)

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

nPat = size(inputArgList, 1);
load(sprintf('%sanalysis/STRFmetrics_%s_%ipat.mat', workingDir, suffix, nPat), 'metrics', 'metricsCols');

idxColSigR = ismember(metricsCols, 'sigR');
idxColR2 = ismember(metricsCols, 'R2');
threshR2 = 0;
idxSig = find(metrics(:, idxColSigR) == 1 & metrics(:, idxColR2) > threshR2 & sum(metrics(:, end-1:end), 2) == 0);

nElecs = length(idxSig);
sigElecsMat = metrics(idxSig, 1:2);

nFreq = 32; nLags = 75;
nFeat = nFreq*nLags;
STRFmat = zeros(nFeat, nElecs);
for idx = 1:nElecs
    load(sprintf('%sencoding/%s_TheWall1_run%i_%s_e%i.mat', workingDir, inputArgList{sigElecsMat(idx, 1), :}, suffix, sigElecsMat(idx, 2)), 'save_coefs');
    STRF = squeeze(mean(save_coefs) ./ std(save_coefs));
    STRFmat(:, idx) = reshape(STRF, 1, nFeat);
end
facTMP = factor(nElecs);
nRows = facTMP(end);
nCols = prod(facTMP(1:end-1));
gridNumbering = reshape(1:nElecs, nRows, nCols);

save(sprintf('%sanalysis/STRFmat_%s_%ipat.mat', workingDir, suffix, nPat), 'STRFmat', 'gridNumbering');