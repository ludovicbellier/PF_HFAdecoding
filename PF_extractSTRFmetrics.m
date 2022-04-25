function [metrics, metricsCols] = PF_extractSTRFmetrics(inputArgList, suffix)

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

nPat = size(inputArgList, 1);
latDict = {'lh', 'rh', 'b'};

metrics = cell(nPat, 1);
idPat = cell(nPat, 1);
idElec = cell(nPat, 1);
latElec = cell(nPat, 1);
typeElec = cell(nPat, 1);
for idxPatient = 1:nPat
    [patientCode, runCode] = inputArgList{idxPatient, :};
    suffixTMP = sprintf('TheWall1_run%i_%s', runCode, suffix);
    [metrics{idxPatient}, metricsCols] = PF_getModelMetrics(patientCode, suffixTMP, 'encoding', workingDir);
    nElecsTMP = size(metrics{idxPatient}, 1);
    idPat{idxPatient} = idxPatient * ones(nElecsTMP, 1);
    idElec{idxPatient} = (1:nElecsTMP)';
    patientInfo = PF_infoPatients(inputArgList{idxPatient});
    latElec{idxPatient} = zeros(nElecsTMP, 1) + find(ismember(latDict, patientInfo.laterality));
    typeElec{idxPatient} = ones(nElecsTMP, 1);
    if ~isempty(patientInfo.ecogStrips)
        typeElec{idxPatient}([patientInfo.ecogStrips{:}]) = 2;
    end
    if ~isempty(patientInfo.ecogDepth)
        typeElec{idxPatient}([patientInfo.ecogDepth{:}]) = 3;
    end
end
metrics = [cat(1, idPat{:}), cat(1, idElec{:}), cat(1, latElec{:}), cat(1, typeElec{:}), cat(1, metrics{:})];
metricsCols = horzcat({'idPat', 'idElec', 'latElec', 'typeElec'}, metricsCols);

fname = sprintf('%sanalysis/STRFmetrics_%s_%ipat.mat', workingDir, suffix, nPat);
save(fname, 'metrics', 'metricsCols');