global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

idxTestSet = 6001:7500;

[reconMatLin, r2actualLin] = PF_plotRecon('supergrid_TheWall1_run1_HFA_params19_f', [-1 0 0]);
[reconMatNonLin, r2actualNonLin] = PF_plotRecon('supergrid_TheWall1_run1_HFA_params96_f', [-1 0 0]);

[H, P, CI, STATS] = ttest(r2actualNonLin, r2actualLin);


%% Fig. 6C
load(sprintf('%s_stimuli/thewall1_stim128.mat', workingDir));
stim128 = stim128(idxTestSet, :)';
tb = linspace(60.01, 75, 1500);
yTicks = [1 12:24:128 128];
yTickLabels = [180 250*2.^(0:4) 7200];

figure('Position', [6 47 1398 1249], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
subplot(311);
imagesc(tb, 1:128, stim128, [0 1]); axis xy; xlim([60 75]);
set(gca, 'XTick', 60:2:75, 'YTick', yTicks, 'YTickLabel', yTickLabels)
ylabel('frequency (Hz)');
xlabel('time (s)');

subplot(312);
imagesc(tb, 1:128, reconMatLin, [0 1]); axis xy; xlim([60 75]);
set(gca, 'XTick', 60:2:75, 'YTick', yTicks, 'YTickLabel', yTickLabels)
ylabel('frequency (Hz)');
xlabel('time (s)');

subplot(313);
imagesc(tb, 1:128, reconMatNonLin, [0 1]); axis xy; xlim([60 75]);
set(gca, 'XTick', 60:2:75, 'YTick', yTicks, 'YTickLabel', yTickLabels)
ylabel('frequency (Hz)');
xlabel('time (s)');


%% Supp. Mat. Audio Files - create wav files
[~, ~, linWav] = PF_plotRecon('supergrid_TheWall1_run1_HFA_params19_f', [-1 0 500]);
audiowrite(sprintf('%saudioRecon/supergrid_TheWall1_run1_params19_bestRsAuto_500iter_DE2.wav', workingDir), linWav/max(abs(linWav))*.95, 16000);

[~, ~, nonLinWav] = PF_plotRecon('supergrid_TheWall1_run1_HFA_params96_f', [-1 0 500]);
audiowrite(sprintf('%saudioRecon/supergrid_TheWall1_run1_params96_bestRsAuto_500iter_DE2.wav', workingDir), nonLinWav/max(abs(nonLinWav))*.95, 16000);

[~, origWav] = PF_aud2wav(stim128', 0, 500);
audiowrite(sprintf('%saudioRecon/original_TheWall1_wav2aud2wav_500iter.wav', workingDir), origWav/max(abs(origWav))*.95, 16000);

[origStim, fs] = audioread(sprintf('%s_stimuli/thewall1.wav', workingDir));
tb = linspace(0, (length(origStim)-1)/fs, length(origStim));
[~, idxBeg] = min(abs(tb-60));
[~, idxEnd] = min(abs(tb-75));
origStim = origStim(idxBeg:idxEnd-1);
audiowrite(sprintf('%saudioRecon/original_TheWall1_60to75s.wav', workingDir), origStim/max(abs(origStim))*.95, fs);