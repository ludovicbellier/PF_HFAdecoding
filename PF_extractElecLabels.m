function PF_extractElecLabels(patientList)

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

labelId = 'fs';
nPat = size(patientList, 1);
labels = cell(nPat, 1);
if strcmp(labelId, 'fs')
    for idxPat = 1:nPat
        fnameLabels = sprintf('%s_anatomy/%s/%s_FsLabels.mat', workingDir, patientList{idxPat}, patientList{idxPat});
        if exist(fnameLabels, 'file')
            load(fnameLabels, 'FsLabels');
            labels{idxPat} = FsLabels;
        else
            load(sprintf('%s_anatomy/%s/%s_labelTable.mat', workingDir, patientList{idxPat}, patientList{idxPat}), 'labelTable_FS', 'labelNames_FS');
            [~, idxLabels] = max(labelTable_FS(:, 2:end), [], 2);
            labels{idxPat} = labelNames_FS(idxLabels);
        end
    end
end
labels = vertcat(labels{:});

save(sprintf('%sanalysis/PFanatLabels_%ipat.mat', workingDir, nPat), 'labels');