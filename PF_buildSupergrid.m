function PF_buildSupergrid

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

load(sprintf('%sanalysis/STRFmetrics_HFA_29pat.mat', workingDir), 'metrics', 'metricsCols');

idxColSigR = ismember(metricsCols, 'sigR');
idxColR2 = ismember(metricsCols, 'R2');
threshR2 = 0;
if ~isnan(threshR2)
    idxSig = find(metrics(:, idxColSigR) == 1 & metrics(:, idxColR2) > threshR2 & sum(metrics(:, end-1:end), 2) == 0);
else
    idxSig = find(metrics(:, idxColSigR) == 1 & sum(metrics(:, end-1:end), 2) == 0);
end
sigElecs = metrics(idxSig, :);
nElecs = size(sigElecs, 1);

patList = unique(sigElecs(:, 1));
nPat = length(patList);
[patientCode, runCode] = inputArgList{patList(1), :};
load(sprintf('%s_preprocessed_ECoG/%s_TheWall1_run%i_preprocessed_HFA.mat', workingDir, patientCode, runCode), 'ecog', 'stim32', 'CF32', 'stim128', 'CF128', 'patientInfo');
L = size(ecog, 1);
supergrid = zeros(L, nElecs);
supergridArtifacts = false(L, nElecs);
counter = 0;
for idxPat = 1:nPat
    [patientCode, runCode] = inputArgList{patList(idxPat), :};
    fprintf('Processing patient %s run %i...\n', patientCode, runCode);
    load(sprintf('%s_preprocessed_ECoG/%s_TheWall1_run%i_preprocessed_HFA.mat', workingDir, patientCode, runCode), 'ecog', 'artifacts');
    idxTMP = sigElecs(sigElecs(:, 1) == patList(idxPat), 2);
    idxTMP2 = counter + (1:length(idxTMP));
    supergrid(:, idxTMP2) = ecog(:, idxTMP);
    supergridArtifacts(:, idxTMP2) = artifacts(:, idxTMP);
    counter = counter + length(idxTMP);
end
ecog = supergrid;
artifacts = supergridArtifacts;

save(sprintf('%s_preprocessed_ECoG/supergrid_TheWall1_run1_preprocessed_HFA.mat', workingDir), 'ecog', 'artifacts', 'stim32', 'CF32', 'stim128', 'CF128', 'patientInfo');
