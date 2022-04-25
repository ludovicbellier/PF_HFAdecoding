global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

% 1A top - stimulus waveform
[y, fs] = audioread(sprintf('%s_stimuli/thewall1.wav', workingDir));
tb = linspace(0, (length(y)/fs)-1, length(y));
figure('Position', [300 600 2000 400], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
plot(tb, y); box off; xlim(tb([1 end]));
print(sprintf('%s_figures/Fig1_stimWaveform.svg', workingDir), '-dsvg');

% 1A bottom - auditory spectrogram
load(sprintf('%s_stimuli/thewall1_stim128.mat', workingDir));
figure('Position', [300 600 2000 500], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
imagesc(stim128'); box off; axis xy; colormap parula;
set(gca, 'CLim', [0 1]);
print(sprintf('%s_figures/Fig1_stimSpecAll.svg', workingDir), '-dsvg');

% 1D top - auditory spectrogram zoom in
xlim([6001 7000]);
print(sprintf('%s_figures/Fig1_stimSpecZoom.svg', workingDir), '-dsvg');

% 1B - X-ray
patientInfo = PF_infoPatients('AMC038');
dcmInfo = dicominfo(patientInfo.Xray{2});
Y = dicomread(dcmInfo);
figure('Position', [6 212 1534 1084], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
imshow(imrotate(imcomplement(Y), 90), [64800 65500]);
print(sprintf('%s_figures/Fig1_XRay.svg', workingDir), '-dsvg');

% 1C - example electrodes
gap = 3;
elecs = 33:36;
gaps = sort(gap:gap:gap*length(elecs), 'descend');
figure('Position', [300 300 2000 800], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
plot(ecog(:, elecs)+gaps, 'color', 'k');
hold on;
plot(ecog(:, 34)+gaps(elecs==34), 'color', 'r');
axis xy; box off;
xlim([1 19072]);
print(sprintf('%s_figures/Fig1_elecExamples.svg', workingDir), '-dsvg');

% 1D bottom - example elec zoom in
load(sprintf('%s_preprocessed_ECoG/AMC038_TheWall1_run1_preprocessed_HFA.mat', workingDir));
figure('Position', [300 600 2000 400], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
plot(ecog(:,34), 'color', 'r'); box off;
xlim([6001 7000]);
print(sprintf('%s_figures/Fig1_elecZoom.svg', workingDir), '-dsvg');

% 1E - STRF
load(sprintf('%s/encoding/AMC038_TheWall1_run1_HFA_e34.mat', workingDir));
STRF = squeeze(mean(saveW2d)) ./ squeeze(std(saveW2d));
figure('Position', [300 300 1000 800], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
imagesc(STRF, [-5 5]); colormap jet; axis xy;
print(sprintf('%s_figures/Fig1_STRFexample.svg', workingDir), '-dsvg');
colorbar;
print(sprintf('%s_figures/Fig1_STRFcolorbar.svg', workingDir), '-dsvg');
