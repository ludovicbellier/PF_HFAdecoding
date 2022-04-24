function PF_preprocessingStim(fname, params, flags)

% Initialize NSLtools
loadload;

% Set default parameters
if nargin < 3
    flags = [1 0]; % save, plot
end
if nargin < 2
    params = [10 10 -2 0];
end
sufParams = ['-' strjoin(arrayfun(@(x) num2str(x), params, 'un', 0), '_')];

% Load the stimulus waveform and resample it for NSLtools
[y, fs] = audioread(fname);
y_16k = resample(y, 16000, fs);

% Compute the 128-fBin auditory spectrogram and characteristic frequencies
stim128 = wav2aud(y_16k, params);
CF128 = 440*2.^((-31:97)/24); % see wav2aud help
CF128 = round(arrayfun(@(i) mean(CF128(i:i+1)), 1:128));

% Plot the auditory spectrogram
if flags(2) == 1
    figure('Position', [6 700 1580 596], 'DefaultAxesFontSize', 5);
    aud_plot(stim128, params);
    title(sprintf('%s%s', fname, sufParams), 'Interpreter', 'none');
end

% Compute the 32-fBin auditory spectrogram and characteristic frequencies
global COCHBA
COCHBA = COCHBA(:, 1:4:end);
stim32 = wav2aud(y_16k, params);
CF32 = 440*2.^((-31:4:97)/24); % see wav2aud help
CF32 = round(arrayfun(@(i) mean(CF32(i:i+1)), 1:32));

% Save both 128- and 32-fBin auditory spectrograms
if max(abs(params - [10 10 -2 0])) == 0
    sufParams = '';
end
if flags(1) == 1
    fname128 = sprintf('%s_stim128%s.mat', fname(1:end-4), sufParams);
    fname32 = sprintf('%s_stim32%s.mat', fname(1:end-4), sufParams);
    save(fname128, 'stim128', 'CF128', 'params');
    save(fname32, 'stim32', 'CF32', 'params');
end