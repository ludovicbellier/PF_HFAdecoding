function PF_preprocessingECoG(inputArgs, workingDir)
% inputArgs = 'AMC062-1';


%% Add FieldTrip to the path and set directories
if exist('ft_defaults.m', 'file') == 0
    ftDir = '/home/knight/lbellier/DataWorkspace/_tools/git/fieldtrip/';
    addpath(ftDir); ft_defaults;
end
if ~exist('workingDir', 'var') || isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end
if ~exist(sprintf('%s_preprocessed_ECoG', workingDir), 'dir')
    mkdir(sprintf('%s_preprocessed_ECoG', workingDir));
end


%% Get input arguments (enables interfacing with SunGridEngine)
inputTMP = split(inputArgs, '-');
patientCode = inputTMP{1};
runCode = str2double(inputTMP{2});


%% Set user-defined parameters
flagFig = 0; % 0 for no figure outputs, 1 for important figures, 2 for all
fsOut = 100; % in Hz
padDur = 10; % padding duration for edge artifacts (each side), in seconds
HPcutoff = 1; % high-pass filter cutoff
LPcutoff = 250; % low-pass filter cutoff
freqBand = [70 150];
sbParams = [20 5 4 1 1]; % bp bandwidth (Hz), step (Hz), order, flagBounded, flagZscore
suffixOut = 'HFA';
trigCode = 1;
stimName = 'thewall1';
screenSetup = 'lab';

stimPrefix = sprintf('TheWall%s_run%i', stimName(end), runCode);
switch screenSetup
    case 'laptop'
        pos1 = [1 581 1530 683];
        pos2 = pos1;
    case 'home'
        pos1 = [1 42 1227 1222];
        pos2 = [1 42 1227 1222];
    otherwise
        pos1 = [1 45 1144 1200];
        pos2 = [1145 45 1144 1200];
end


%% Get patient information
fprintf('Processing patient %s...\n', patientCode);
patientInfo = PF_infoPatients(patientCode);
nElecs = length(patientInfo.ecogElec{1});
refElec = patientInfo.refElec{runCode};
noisyElecs = patientInfo.noisyElecs{runCode};
epilepticElecs = patientInfo.epilepticElecs{runCode};
filename = patientInfo.filename{runCode};
CARsubsets = patientInfo.CARsubsets{runCode};
goodElecs = setdiff(1:nElecs, [refElec noisyElecs epilepticElecs]);
goodElecsDisp = setdiff(1:nElecs, noisyElecs);


%% Load ECoG data
cfg = [];
cfg.channel = patientInfo.ecogElec{1};
cfg.dataset = filename;
dataRaw = ft_preprocessing(cfg);
fs = dataRaw.fsample;
[~, triggers] = load_bcidat(cfg.dataset);


%% Visualize whole raw data
if flagFig > 0
    cfg = [];
    cfg.viewmode = 'vertical';
    cfg.preproc.demean = 'yes';
    cfg.blocksize = 10;
    cfg.channel = goodElecsDisp;
    cfg.position = pos1;
    cfg.ylim = [-1 1]*400;
    ft_databrowser(cfg, dataRaw);
end


%% Get stimulus length
load(sprintf('%s_stimuli/%s_stim32.mat', workingDir, stimName), 'stim32');
L = (size(stim32, 1) / fsOut) * fs;


%% Get triggers
stimBeg = find(triggers.StimulusCode==trigCode, 1, 'first');
if stimBeg == 1 % avoids spurious start
    stimBeg = find(triggers.StimulusCode(5000:end)==trigCode, 1, 'first') + 4999;
end
stimEnd = stimBeg + L - 1;


%% Plot trigger channel
if flagFig == 2
    figure('Position', [110 510 725 570], 'DefaultAxesFontSize', 5); hold on;
    plot(triggers.StimulusCode);
    plot(stimBeg, trigCode, '*r');
    plot(stimEnd, trigCode, '*r');
    xlabel('time (samples)');
    ylabel('stimulus code');
    xlim([1 length(triggers.StimulusCode)]);
    ylim([-.1 max(double(triggers.StimulusCode))+.1]);
    title(sprintf('%s %s %i - stimBeg=%i - stimEnd=%i', patientCode, stimName, runCode, stimBeg, stimEnd), 'Interpreter', 'none');
    input('Press any key to proceed...\n');
    close all
end


%% Cut data with data padding
padding = padDur * fs;
if stimEnd + padding > size(dataRaw.trial{1}, 2) % mirror data if insufficient length for padding
    dataRaw.time{1}(end+1:end+padding) = dataRaw.time{1}(end) + dataRaw.time{1}(2:padding+1);
    dataRaw.trial{1}(:, end+1:end+padding) = dataRaw.trial{1}(:, end:-1:end-padding+1);
    dataRaw.sampleinfo = [1 length(dataRaw.time{1})];
end
cfg = [];
cfg.begsample = stimBeg - padding;
cfg.endsample = stimEnd + padding;
data = ft_redefinetrial(cfg, dataRaw);
Lpad = size(data.trial{1}, 2);


%% Visualize raw data
if flagFig > 0
    cfg = [];
    cfg.viewmode = 'vertical';
    cfg.preproc.demean = 'yes';
    cfg.blocksize = 10;
    cfg.channel = goodElecsDisp;
    cfg.position = pos1;
    cfg.ylim = [-1 1]*400;
    ft_databrowser(cfg, data);
end


%% Choose CAR strategy (on all electrodes vs. on subsets)
if iscell(CARsubsets)
    elecSubsets = CARsubsets;
else
    fnameCARstrat = sprintf('%s_preprocessed_ECoG/%s_%s_%s_covarianceMatrices.png', workingDir, patientCode, stimPrefix, suffixOut);
    elecSubsets = PF_chooseCARstrategy(data, noisyElecs, goodElecs, freqBand, fnameCARstrat);
end


%% Filter out power line noise
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = (60 * (1:round((LPcutoff+60)/60)) + [-1; 1])';
cfg.hpfilter = 'yes';
cfg.hpfreq = HPcutoff;
cfg.lpfilter = 'yes';
cfg.lpfreq = LPcutoff;
data = ft_preprocessing(cfg, data);


%% Visualize clean data
if flagFig > 0
    cfg = [];
    cfg.viewmode = 'vertical';
    cfg.preproc.demean = 'yes';
    cfg.blocksize = 10;
    cfg.channel = goodElecsDisp;
    cfg.position = pos2;
    cfg.ylim = [-1 1]*200;
    ft_databrowser(cfg, data);
end


%% Filter HFA subbands, process CAR and compute Hilbert transform
if sbParams(4) == 0
    subbands = (freqBand(1):sbParams(2):freqBand(2))'+[-1 1]*(sbParams(1)/2);
else
    subbands = (freqBand(1):sbParams(2):freqBand(2)-sbParams(1))'+[0 1]*sbParams(1);
end
nSubbands = size(subbands, 1);
dataSave = zeros(nElecs, Lpad, nSubbands);
for idxSubband = 1:nSubbands
    fprintf('Processing subband #%i/%i...\n', idxSubband, nSubbands);
    
    % 1. bandpass filter / get subband
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [subbands(idxSubband, 1) subbands(idxSubband, 2)];
    cfg.bpfiltord = sbParams(3);
    dataTMP = ft_preprocessing(cfg, data);
    
    % 2. CAR / denoise data
    cfg = [];
    cfg.reref = 'yes';
    cfg.refmethod = 'median';
    if isempty(elecSubsets)
        cfg.refchannel = goodElecs;
        dataTMP = ft_preprocessing(cfg, dataTMP);
    else
        dataCAR = dataTMP.trial{1};
        for j = 1:length(elecSubsets)
            cfg.refchannel = setdiff(elecSubsets{j}, [refElec, noisyElecs, epilepticElecs]);
            dataCARpartial = ft_preprocessing(cfg, dataTMP);
            dataCAR(elecSubsets{j}, :) = dataCARpartial.trial{1}(elecSubsets{j}, :);
        end
        dataTMP.trial{1} = dataCAR;
    end
    
    % 3. Hilbert transform / get envelope
    cfg = [];
    cfg.hilbert = 'abs';
    dataTMP = ft_preprocessing(cfg, dataTMP);
    
    % save subband data
    dataSave(:, :, idxSubband) = dataTMP.trial{1};
end

% normalize subbands before averaging
if sbParams(5) == 1
    if nSubbands > 1
        for idxElec = 1:nElecs
            dataTMP.trial{1}(idxElec, :) = mean(robustScaler(squeeze(dataSave(idxElec, :, :))), 2);
        end
    else
        dataTMP.trial{1} = robustScaler(dataSave, 2);
    end
else
    if nSubbands > 1
        dataTMP.trial{1} = mean(dataSave, 3);
    else
        dataTMP.trial{1} = dataSave;
    end
end


%% Visualize preprocessed data
if flagFig > 0
    cfg = [];
    cfg.viewmode = 'vertical';
    cfg.preproc.demean = 'yes';
    cfg.blocksize = 10;
    cfg.channel = goodElecsDisp;
    cfg.position = pos1;
    cfg.ylim = [-1 1]*1;
    ft_databrowser(cfg, dataTMP);
end


%% Remove padding
cfg = [];
cfg.begsample = padding + 1;
cfg.endsample = padding + L;
data = ft_redefinetrial(cfg, dataTMP);


%% Downsample data
cfg = [];
cfg.resamplefs = fsOut;
data = ft_resampledata(cfg, data);


%% Extract data from its FieldTrip structure
ecog = data.trial{1}';
tb = data.time{1};


%% Save preprocessed data
patientInfo.time{1} = tb;
patientInfo.CARsubsets{runCode} = elecSubsets;

load(sprintf('%s_stimuli/%s_stim32.mat', workingDir, stimName), 'stim32', 'CF32');
load(sprintf('%s_stimuli/%s_stim128.mat', workingDir, stimName), 'stim128', 'CF128');

save(sprintf('%s_preprocessed_ECoG/%s_%s_preprocessed_%s.mat', workingDir, patientCode, stimPrefix, suffixOut), 'stim32', 'stim128', 'CF32', 'CF128', 'ecog', 'patientInfo');