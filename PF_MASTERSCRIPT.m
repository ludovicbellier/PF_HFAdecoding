%{
This script performs all preprocessing and analysis steps underlying our
Pink Floyd | HFA decoding study.

I. Preprocessing:
- preprocess the song stimulus (waveform -> auditory spectrogram);
- for each patient:
    - preprocess ECoG data (estimation of HFA);
    - detect artifacts (noisy time samples);
    - process anatomical data (MNI electrode coordinates and atlas labels);

II. Encoding models:
- launch all encoding models (STRFs);
- collect STRF metrics and identify significant electrodes
- Fig. 1 - plot example data
- Fig. 2 - analyze STRF prediction accuracies
- Fig. 3 - perform an ICA on STRF coefficients
- Fig. 4 - perform a sliding correlation between ICA components and
           the song spectrogram

III. Decoding models:
- launch linear decoding models for all significant electrodes, and after
  ablating anatomical and functional sets of electrodes
- launch linear decoding models by bootstrapping the number of electrodes
  used as predictors and the dataset duration
- launch linear and nonlinear reconstruction of the song spectrogram from
  all significant electrodes
- Fig. 5 - analyze decoding accuracies in the ablation analysis
- Fig. 6 - analyze output of the bootstrap analysis and of the linear and
           nonlinear song reconstruction
%}


%% Initialization
% Define user-specific inputs
rootDir = '/home/knight/lbellier/DataWorkspace/';
mainToolsDir = [rootDir '_tools/git/PF_HFAdecoding/'];
IPADir = [rootDir '_tools/git/IPA/'];
NSLDir = [rootDir '_tools/nsltools/'];
FTDir = [rootDir '_tools/git/fieldtrip/'];

global workingDir;
workingDir = [rootDir '_projects/PinkFloyd/'];
global dataDir;
dataDir = '/home/knight/ecog/DATA_FOLDER/Albany/';

inputArgList = {'AMC006', 1; 'AMC007', 1; 'AMC008', 1; 'AMC009', 1; 'AMC010', 1;...
    'AMC011', 1; 'AMC012', 2; 'AMC013', 1; 'AMC014', 1; 'AMC015', 1;...
    'AMC017', 1; 'AMC018', 2; 'AMC019', 1; 'AMC022', 1; 'AMC026', 1;...
    'AMC027', 1; 'AMC028', 1; 'AMC029', 3; 'AMC031', 1; 'AMC032', 1;...
    'AMC033', 1; 'AMC037', 1; 'AMC038', 1; 'AMC039', 2; 'AMC040', 2;...
    'AMC041', 2; 'AMC044', 1; 'AMC045', 1; 'AMC062', 1};

suffix = 'HFA';

% Add toolboxes and functions to the path
pathDirs = {mainToolsDir, IPADir, NSLDir, FTDir};
addpath(pathDirs{:});


%% I. Preprocessing
% Preprocess the song stimulus (waveform -> auditory spectrogram)
fnameStim = [workingDir '_stimuli/thewall1.wav'];
if ~exist(fnameStim, 'file')
    PF_preprocessingStim(fnameStim);
end

% for each patient
nPat = length(inputArgList);
for idxPat = 1:nPat
    [patientCode, runCode] = inputArgList{idxPat, :};
    
    % first, fill up information in PF_infoPatients.m
    % edit PF_infoPatients.m

    % Preprocess ECoG data (estimation of HFA)
    PF_preprocessingECoG(sprintf('%s-%i', patientCode, runCode));
    % PF_SGElauncher_preprocessing(patientCode, runCode);
    
    % Detect artifacts (noisy time samples)
    PF_preprocessingArtifacts(patientCode, runCode, suffix, 'manual');

    % Process anatomical data (MNI electrode coordinates and atlas labels)
    PF_processingAnatomy(patientCode);
end


%% Encoding analysis
% Launch model fittings
% suffixList = {'LFC', '1-4Hz', '4-8Hz', '8-13Hz', '13-30Hz', '30-70Hz', '70-150Hz'};
suffixList = {'HFB'};
nSuf = length(suffixList);
for idxPat = 1:nPat
    [patientCode, runCode] = inputArgList{idxPat, :};
    for idxSuf = 1:nSuf
        PF_SGElauncher_modeling(patientCode, runCode, 'encoding', suffixList{idxSuf});%, [], 1);
    end
end

% Check progress
fprintf('STRF jobs progress\n-----------------\n');
for idxPat = 1:nPat
    patientCode = inputArgList{idxPat, 1};
    infoPat = PF_infoPatients(patientCode);
    nElecs = length(infoPat.ecogElec{1});
    nSTRFs = length(dir(sprintf('encoding/%s*%s*mat', patientCode, suffix)));
    fprintf('%2i - %s - %3i/%3i - %3i%%\n', idxPat, patientCode, nSTRFs, nElecs, round(100*nSTRFs/nElecs));
end

% PF_buildElecInfoMatrix;

% Output plots
idxElec = 70;
filename = sprintf('encoding/%s_TheWall1_run%i_%s_e%i.mat', patientCode, runCode, suffix, idxElec);
PF_plotSTRF(filename);

figOut = 0;
for idxPat = 1:nPat
    [patientCode, runCode] = inputArgList{idxPat, :};
    superSuffix = sprintf('TheWall1_run%i_%s', runCode, suffix);
    PF_plotSTRFs(patientCode, superSuffix, 'r', figOut);
    PF_plotSTRFmetrics(patientCode, superSuffix, 'r', figOut);
    PF_plotValOnAnat({patientCode, runCode}, suffix, 'sigR', 'patient', 2, figOut);
    input('');
    close all;
end

idxElec = 56;
for idxSuf = 1:length(suffixList)
    filename = sprintf('encoding/%s_TheWall1_run%i_%s_e%i.mat', patientCode, runCode, suffixList{idxSuf}, idxElec);
    PF_plotSTRF(filename);
end


%% Decoding analysis
PF_SGElauncher_modeling(patientCode, runCode, 'decoding', 'HFB');

% Output plots
superSuffix = sprintf('TheWall1_run%i_%s', runCode, suffix);
PF_plotDecodingMetrics(patientCode, superSuffix, 'decoding');

% Ablations
PF_ablationAnalysis


%% Recon analysis
PF_SGElauncher_modeling('supergrid', 1, 'recon', 'HFB', setdiff(1:128, 38), 1); % params1 - linear regression
PF_SGElauncher_modeling('supergridSTGwoSus', 1, 'recon', 'HFB', 1:107, 2); % params2 - linear regression wo early stopping
PF_SGElauncher_modeling('supergridSTGwoSus', 1, 'recon', 'HFB', 1:107*2, 3); % params3 - ridge regression

patientCode = 'supergridSTGwoSus'; runCode = 1; modelType = 'recon'; suffix = 'HFB'; paramCode = 1;
patientCode = 'supergrid'; runCode = 1; modelType = 'recon'; suffix = 'HFB'; paramCode = 7;
PF_SGElauncher_modeling(patientCode, runCode, modelType, suffix, 1:107, paramCode);

% params1: rMLRwES; params2: rMLRwoES; params3: ridge(woES);
% params4: MLP; params5: MLPwoES; params6: MLPw(HLS64,64)(a10); params7: MLPw(HLS64,64)(a1)

% [S, S80best] = PF_recon2aud('supergridSTGwoSus', 'TheWall1_run1_HFB_params1', 'recon', 1);
paramsRecon = [50 0 1 0]; % threshPrctile, threshMetrics, flagFig, forceRecon
for paramCode = 1:7
    recon = PF_plotRecon(sprintf('%s_TheWall1_run%i_%s_params%i_f', patientCode, runCode, suffix, paramCode), paramsRecon);
    [wav0, wav, err] = PF_aud2wav(recon', 0, 250); % recon' must be L*nFreq
    save(sprintf('audioRecon/%s_TheWall1_run%i_%s_params%i_recon.mat', patientCode, runCode, suffix, paramCode), 'wav0', 'wav', 'err')
    pause(1);
    audiowrite(sprintf('audioRecon/%s_TheWall1_run%i_%s_params%i_recon.wav', patientCode, runCode, suffix, paramCode), wav./max(abs(wav)), 16000);
end

paramsRecon = [50 1 0 1]; % threshPrctile, threshMetrics, flagLogCorrect, flagFig
recon = PF_plotRecon(sprintf('%s_TheWall1_run%i_%s_params%i_f', patientCode, runCode, suffix, paramCode), paramsRecon);
test1 = smoothdata(recon, 'gaussian', 5);
test2 = conv2(recon, (1/9)*ones(3), 'same');
[~, wav] = PF_aud2wav(recon', 0, 5); % recon' must be L*nFreq
soundsc(wav, 16000);


fname = 'recon/supergrid_TheWall1_run1_HFB_params4_f2.mat';
PF_plotReconMLPcv(fname, 2);

nParamCodes = 7;
metrics = cell(nParamCodes, 1);
colTMP = lines(nParamCodes);
for paramCode = 1:nParamCodes
    [~, metrics{paramCode}] = PF_plotRecon(sprintf('%s_TheWall1_run%i_%s_params%i_f', patientCode, runCode, suffix, paramCode), 0);
    plot_meanWithShadedSEM(metrics{paramCode}(:, 3), metrics{paramCode}(:, 4), 1:128, [], colTMP(paramCode, :));
end
legendStr = sprintf('''%i'','''','''','''',', 1:nParamCodes);
eval(sprintf('legend(%s);', legendStr(1:end-1)));

% recon from regular, variable test set
load decoding/supergrid_TheWall1_run1_HFB_params50_f3.mat
L = max(cellfun(@max, save_split_indices(:, 3)));
N = length(save_y_pred);
sparseRecon = nan(L+1, N);
for idx = 1:N
    sparseRecon(save_split_indices{idx, 3}+1, idx) = save_y_pred{idx};
end
figure('Position', [54 385 2320 825], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
subplot(3,1,1:2); imagesc(sparseRecon');
subplot(313); plot(mean(sparseRecon, 2, 'omitnan'));

% patientCode = 'AMC062'; runCode = 1; suffix = '70-150Hz_params5';
% load(sprintf('analysis/%s_sigElecs.mat', patientCode), 'sigElecs');
% PF_SGElauncher_modeling(patientCode, runCode, 'recon', suffix, sigElecs, []);
% PF_plotDecodingModel('recon/supergrid_TheWall1_run1_HFB_f13.mat');

% superSaveR = zeros(128, 4);
% superSaveR2 = zeros(128, 4);
% for idx = 1:128
%     disp(idx);
%     [superSaveR(idx, 1), superSaveR2(idx, 1)] = PF_plotMLP(sprintf('recon/AMC062_TheWall1_run1_70-150Hz_f%i.mat', idx), 0);
%     [superSaveR(idx, 2), superSaveR2(idx, 2)] = PF_plotMLP(sprintf('recon/AMC062_TheWall1_run1_70-150Hz_20HzBins_f%i.mat', idx), 0);
%     [superSaveR(idx, 3), superSaveR2(idx, 3)] = PF_plotMLP(sprintf('recon/AMC062_TheWall1_run1_70-150Hz_10HzBins_f%i.mat', idx), 0);
%     [superSaveR(idx, 4), superSaveR2(idx, 4)] = PF_plotMLP(sprintf('recon/AMC062_TheWall1_run1_70-150Hz_5HzBins_f%i.mat', idx), 0);
% end
% figure;
% subplot(2,3,1:2); plot(superSaveR); xlim([1 128]);
% subplot(2,3,3); bar(median(superSaveR));
% subplot(2,3,4:5); plot(superSaveR2); xlim([1 128]);
% subplot(2,3,6); bar(median(superSaveR2));

% suffixList = {'70-150Hz_params5', '70-150Hz_params6', '70-150Hz_params7', '70-150Hz_params8'};
% for idxSuf = 1:4
%     suffix = suffixList{idxSuf};
%     disp(suffix);
%     [S, S80best] = PF_recon2aud('AMC062', 1, suffix, 1);
%     [~, y] = PF_aud2wav(S, 1, 100);
%     pause(1);
%     audiowrite(sprintf('recon/%s_TheWall1_run%i_%s_recon.wav', patientCode, runCode, suffix), y./max(abs(y)), 16000);
%     [~, y80best] = PF_aud2wav(S80best, 1, 100);
%     pause(1);
%     audiowrite(sprintf('recon/%s_TheWall1_run%i_%s_recon80best.wav', patientCode, runCode, suffix), y80best./max(abs(y80best)), 16000);
% end

% new 5/2020
% suffix2List = {'-fac-2_params3', '-fac-2_params4'}; modelType = 'recon';
% for idx = 1:length(suffix2List)
%     suffix2 = suffix2List{idx};
%     superSuffix = sprintf('TheWall1_run%i_%s%s', runCode, suffix, suffix2);
%     PF_plotDecodingMetrics(patientCode, superSuffix, modelType);
%     [recon, reconXpBest] = PF_recon2aud(patientCode, superSuffix, modelType, 1);
%     [~, ymin, err] = PF_aud2wav(recon, 1, 100);
%     [~, yminXpBest, errXpBest] = PF_aud2wav(reconXpBest, 1, 100);
%     audiowrite(sprintf('audioRecon/%s_%s_%s_100iters.wav', patientCode, superSuffix, modelType), ymin/max(abs(ymin)), 16000);
%     audiowrite(sprintf('audioRecon/%s_%s_%s_95pBest_100iters.wav', patientCode, superSuffix, modelType), yminXpBest/max(abs(yminXpBest)), 16000);
% end
% idxFbin = 17;
% for idx = 1:length(suffix2List)
%     suffix2 = suffix2List{idx};
%     filename = sprintf('%s/%s_TheWall1_run%i_%s%s_f%i.mat', modelType, patientCode, runCode, suffix, suffix2, idxFbin);
%     PF_plotDecodingModel(filename);
% end



%% Analysis
suffix = '70-150Hz';
inputArgList = {'AMC006', 1; 'AMC007', 1; 'AMC008', 1; 'AMC009', 1; 'AMC010', 1;...
    'AMC011', 1; 'AMC012', 2; 'AMC013', 1; 'AMC014', 1; 'AMC015', 1;...
    'AMC017', 1; 'AMC018', 2; 'AMC019', 1; 'AMC022', 1; 'AMC026', 1;...
    'AMC027', 1; 'AMC028', 1; 'AMC029', 3; 'AMC031', 1; 'AMC032', 1;...
    'AMC033', 1; 'AMC037', 1; 'AMC038', 1; 'AMC039', 2; 'AMC040', 2;...
    'AMC041', 2; 'AMC044', 1; 'AMC045', 1; 'AMC062', 1};
PF_plotValOnAnat(inputArgList, suffix, 'sigR', 'MNI', 2);

% covariance between metrics
[metrics, metricsList] = PF_getModelMetrics(patientCode, ['TheWall1_run1_' suffix], 'encoding');
idxMetrics = [1 2 5:7 9:11];
figure('Position', [100 500 820 650]); imagesc(cov(zscore(metrics(:, idxMetrics)))); colorbar; colormap(redblue);
set(gca, 'CLim', [-1 1], 'YTickLabel', metricsList(idxMetrics), 'XTickLabel', metricsList(idxMetrics), 'XTickLabelRotation', 90);

% STRF measurements
patientCode = 'AMC062'; runCode = 1; suffix = 'HFB';
load analysis/AMC062_sigElecs.mat
nElecs = length(sigElecs);
idxMeas = [1:16 18 20];
valMat = zeros(nElecs, length(idxMeas));
for idxElec = 1:nElecs
    measTMP = PF_getSTRFmeasurements(patientCode, runCode, suffix, sigElecs(idxElec), 1.96);
    valMat(idxElec, :) = cell2mat(measTMP(idxMeas, 2));
    if idxElec == 1
        valNames = measTMP(idxMeas, 1)';
    end
end
% valToPlot = 'entrSigLF';
% valIdx = find(strcmp(valNames, valToPlot));
valIdx = 0;
valIdx = valIdx+1; fprintf('plotting %s...\n', valNames{valIdx});
values = zeros(250, 1);
values(sigElecs) = valMat(:, valIdx);
values(sum(metrics(:, 16:17), 2) > 0) = nan; % don't plot noisy/epi elecs
% values(76)=nan;%%TMP
PF_plotValOnAnat({patientCode, runCode}, suffix, values, 'patient', 2, 0, 'rh');%, [nan 500]);
set(gcf, 'Position', [1 1 1126 970]);

load _anatomy/AMC062/AMC062_gridNumbering.mat
figure('Position', [488 264 1204 416]);
imagesc(values(gridNumbering{1})); colorbar;

for idxVal = 1:length(valNames)
    values = zeros(250, 1);
    values(sigElecs) = valMat(:, idxVal);
    values(sum(metrics(:, 13:14), 2) > 0) = nan; % don't plot noisy/epi elecs
%     margin = values(values>0);
%     margin = [prctile(margin, 20) prctile(margin, 80)];
    subplot(5,4,idxVal); imagesc(values(gridNumbering{1})); colorbar;
    title([num2str(idxVal) ' - ' valNames{idxVal}]);
end

% rhythmic STRFs - run exp_PF_PCAthenICAonSTRFs.m first until STRFmat
% rhythmicSTRFs = [6 8; 8 21; 10 39; 11 70; 12 85; 14 20; 14 88; 14 91;
%     14 94; 15 50; 17 52; 18 86; 19 40; 26 31; 26 73; 26 69; 26 66; 26 54;
%     26 59; 26 64; 27 32; 27 17; 29 57; 29 58; 31 9; 31 83; 37 4; 37 68;
%     37 74; 37 69; 37 86; 37 111; 38 24; 38 23; 38 26; 38 35; 38 34; 40 34;
%     41 28; 44 15; 45 40; 45 41; 45 30; 62 39; 62 59; 62 57; 62 184; 62 190;
%     62 191; 62 192; 62 187; 62 197; 62 79; 62 78; 62 77; 62 195; 62 193; 62 206];
rhythmicSTRFs = [6 8; 8 21; 11 70; 12 85; 14 20; 14 88; 14 91;
    14 94; 26 73; 26 69; 26 66; 
    26 64; 29 57; 31 9; 31 83; 37 4;
    37 86; 37 111; 38 24; 38 26; 38 35; 38 34; 
    41 28; 44 15; 45 41; 45 30; 62 59; 62 57; 62 184; 62 190;
    62 191; 62 192; 62 197; 62 78; 62 77; 62 195; 62 193; 62 206];
idxPatients = cellfun(@(x) str2double(x(4:6)), inputArgList(:, 1));
[~, rhythmicSTRFs(:, 1)] = ismember(rhythmicSTRFs(:, 1), idxPatients);
[~, idx] = ismember(rhythmicSTRFs, sigElecsMat, 'rows');
idx = idx(idx~=0);
N = length(idx);

for i = 1:N
    subplot(6, 10, i);
    STRFtmp = reshape(STRFmat(:, idx(i)), 32, 75);
    imagesc(STRFtmp, [-5 5]); axis xy;
    title(sprintf('%s e%i', inputArgList{sigElecsMat(idx(i), 1), 1}, sigElecsMat(idx(i), 2)));
end


%% FIGURES
PF_extractSTRFmetrics(inputArgList, suffix, workingDir);
PF_extractElecLabels(inputArgList(:, 1), suffix, workingDir);
PF_buildSTRFmat(inputArgList, suffix, workingDir);

PF_fig1_plotStimAndHFB

PF_fig2_getHistAnatSigElecs
PF_fig2explo_STRFmeasurements

PF_fig3_STRFsICArhythm
exp_PF_PCAthenICAonSTRFs
exp_getTemporalModulations

PF_fig4_corrICAcompWithAudSpctrgrm

PF_fig5_ablation

PF_SGElauncher_bootstrap
PF_decodingBootstrap