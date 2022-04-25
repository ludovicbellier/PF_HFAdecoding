function [metrics, metricsCols] = PF_getModelMetrics(patientCode, suffix, modelType)

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

pVal = .05;
threshCI = 1 - 2 * pVal;
threshR2 = 0;

% get elec info
if strcmp(modelType, 'encoding')
    patientInfo = PF_infoPatients(patientCode);
    nChans = length(patientInfo.ecogElec{1});
    runCode = str2double(suffix(13));
    refElec = patientInfo.refElec{runCode};
    noisyElecs = patientInfo.noisyElecs{runCode};
    epilepticElecs = patientInfo.epilepticElecs{runCode};
    elecInfo = zeros(nChans, 3);
    elecInfo(refElec, 1) = 1;
    elecInfo(noisyElecs, 2) = 1;
    elecInfo(epilepticElecs, 3) = 1;
    prefixChan = 'e';
else
    prefixChan = 'f';
    load(sprintf('%s%s/%s_%s_%s2.mat', workingDir, modelType, patientCode, suffix, prefixChan), 'params');
    nChans = str2double(params.stim_type(5:end));
end

% get model metrics (r, R2, sig, sigR2, iter, time, nanCount, isNan, maxCount, LR)
metricsCols = {'r', 'R2', 'sigR', 'sigR2', 'iter', 'time', 'nanCount', 'isNan', 'maxCount', 'LR'};
nMetrics = length(metricsCols);
metrics = zeros(nChans, nMetrics);
for idx = 1:nChans
    try
        fnameTMP = sprintf('%s%s/%s_%s_%s%i.mat', workingDir, modelType, patientCode, suffix, prefixChan, idx);
        assert(exist(fnameTMP, 'file') == 2);
        load(fnameTMP, 'save_r', 'save_r2', 'save_iter', 'totalTime', 'nan_counter', 'max_iter_counter', 'params');
        
        idxNan = find(isnan(save_r));
        if ~isempty(idxNan)
            idxRsOk = 1:idxNan-1;
        else
            idxRsOk = 1:length(save_r);
        end
        metrics(idx, 1) = mean(save_r(idxRsOk));
        metrics(idx, 2) = mean(save_r2(idxRsOk));
        metrics(idx, 3) = double(checkCorrCoefSignificance(save_r(idxRsOk), threshCI));
        metrics(idx, 4) = metrics(idx, 2) > threshR2;
        metrics(idx, 5) = mean(save_iter(idxRsOk));
        metrics(idx, 6) = totalTime;
        metrics(idx, 7) = nan_counter;
        metrics(idx, 8) = double(~isempty(idxNan));
        metrics(idx, 9) = max_iter_counter;
        metrics(idx, 10) = params.learning_rate;
    catch
        fprintf('/!\\ WARNING /!\\ "%s_%s_%s%i.mat" not found\n', patientCode, suffix, prefixChan, idx);
        metrics(idx, :) = ones(1, nMetrics) .* -1;
    end
    clear save_r
end

if strcmp(modelType, 'encoding')
    metrics = [metrics elecInfo];
    metricsCols = [metricsCols 'isRef' 'isNoisy' 'isEpi'];
end