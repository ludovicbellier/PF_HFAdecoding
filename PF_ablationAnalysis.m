function PF_ablationAnalysis

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

load(sprintf('%sanalysis/PFanatLabels_29pat.mat', workingDir), 'labels');
load(sprintf('%sanalysis/STRFmetrics_HFA_29pat.mat', workingDir), 'metrics');
load(sprintf('%sanalysis/STRFfuncSets.mat', workingDir), 'idxFuncSets');
idxSig = metrics(:, 7) == 1 & metrics(:, 6) > 0 & sum(metrics(:, 17:18), 2) == 0;
metrics = metrics(idxSig, :);
labels = labels(idxSig);
N = length(labels);

idxLeft = find(metrics(:, 3) == 1);
idxRight = find(metrics(:, 3) == 2);
idxSTG = find(contains(labels, 'superiortemporal'));
idxSMC = find(contains(labels, 'central'));
idxIFG = find(contains(labels, 'pars'));
idxOther = setdiff(1:N, [idxSTG; idxSMC; idxIFG])';
idxOnset = find(idxFuncSets(:, 1));
idxSustained = find(idxFuncSets(:, 2));
idxLateOnset = find(idxFuncSets(:, 3));
idxRhythmic = find(idxFuncSets(:, 4));
idxList = {[], idxLeft, idxRight,...
    idxSTG, intersect(idxSTG, idxLeft), intersect(idxSTG, idxRight),...
    idxSMC, intersect(idxSMC, idxLeft), intersect(idxSMC, idxRight),...
    idxIFG, intersect(idxIFG, idxLeft), intersect(idxIFG, idxRight), idxOther,...
    idxOnset, intersect(idxOnset, idxLeft), intersect(idxOnset, idxRight),...
    idxSustained, intersect(idxSustained, idxLeft), intersect(idxSustained, idxRight),...
    idxLateOnset, intersect(idxLateOnset, idxLeft), intersect(idxLateOnset, idxRight),...
    idxRhythmic, intersect(idxRhythmic, idxLeft), intersect(idxRhythmic, idxRight)};

% subtractive models (ablation per se - all but x)
suffList = {'', 'NoL', 'NoR',...
            'NoSTG', 'NoLSTG', 'NoRSTG',...
            'NoSMC', 'NoLSMC', 'NoRSMC',...
            'NoIFG', 'NoLIFG', 'NoRIFG', 'NoOther',...
            'NoOnset', 'NoLOnset', 'NoROnset',...
            'NoSustained', 'NoLSustained', 'NoRSustained',...
            'NoLateOnset', 'NoLLateOnset', 'NoRLateOnset',...
            'NoRhythmic', 'NoLRhythmic', 'NoRRhythmic'};
nSuf = length(suffList);

for idx = 1:nSuf
    load(sprintf('%s_preprocessed_ECoG/supergrid_TheWall1_run1_preprocessed_HFA.mat', workingDir), 'ecog', 'artifacts', 'stim32', 'CF32', 'stim128', 'CF128', 'patientInfo');
    ecog = ecog(:, setdiff(1:N, idxList{idx}));
    artifacts = artifacts(:, setdiff(1:N, idxList{idx}));
    save(sprintf('%s_preprocessed_ECoG/supergrid%s_TheWall1_run1_preprocessed_HFA.mat', workingDir, suffList{idx}), 'ecog', 'artifacts', 'stim32', 'CF32', 'stim128', 'CF128', 'patientInfo');

    PF_SGElauncher_modeling(sprintf('supergrid%s', suffList{idx}), 1, 'decoding', 'HFA', 1:size(ecog, 2));
end
