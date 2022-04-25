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
- Fig. 5 - analyze decoding accuracies in the ablation analysis
- launch linear decoding models by bootstrapping the number of electrodes
  used as predictors and the dataset duration
- launch linear and nonlinear reconstruction of the song spectrogram from
  all significant electrodes
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
    PF_preprocessingECoG(sprintf('%s-%i', patientCode, runCode), workingDir);
    % PF_SGElauncher_preprocessing(patientCode, runCode);
    
    % Detect artifacts (noisy time samples)
    PF_preprocessingArtifacts(patientCode, runCode, suffix, 'manual');

    % Process anatomical data (MNI electrode coordinates and atlas labels)
    PF_processingAnatomy(patientCode);
end

% Once all anatomical recons are done, gather anatomical labels
PF_extractElecLabels(inputArgList(:, 1));


%% II. Encoding models
% Launch all encoding models
for idxPat = 1:nPat
    [patientCode, runCode] = inputArgList{idxPat, :};
    PF_SGElauncher_modeling(patientCode, runCode, 'encoding', suffix);
end

% Once all models have run, gather results (metrics + STRFs)
PF_extractSTRFmetrics(inputArgList, suffix);
PF_buildSTRFmat(inputArgList, suffix);

% Get figure elements for Fig. 1, 2 and 3
PF_fig1_plotStimAndHFA;
PF_fig2_getHistAnatSigElecs;
PF_fig3_STRFexamples;

% Perform ICA on STRF coefficients
PF_fig3_ICAonSTRFs;

% Get temporal modulation indices
PF_fig3_temporalModulations;

% Perform sliding correlation between ICA components and song spectrogram
PF_fig4_corrICAcompWithAudSpctrgrm;


%% III. Decoding models
% Build supergrid (all 347 electrodes with a sig. STRF)
PF_buildSupergrid;

% Perform the ablation analysis
PF_ablationAnalysis;

% Analyze results of the ablation analysis
PF_fig5_ablation;

% Launch bootstrap analysis of nb of elecs and dataset duration
PF_SGElauncher_bootstrap;
PF_SGElauncher_bootstrapLength;

% Analyze results of the bootstrap analyses
PF_fig6_decodingBootstrap;
PF_fig6_decodingBootstrapLength;

% Linear vs nonlinear reconstruction
PF_SGElauncher_modeling('supergrid', 1, 'recon', 'HFA', 1:128, 19); % params19 - linear regression
PF_SGElauncher_modeling('supergrid', 1, 'recon', 'HFA', 1:128, 96); % params96 - MLP / nonlinear regression

% Plot recons and create wav files
PF_fig6_recon;
