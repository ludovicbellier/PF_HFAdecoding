function PF_preprocessingArtifacts(patientCode, runCode, suffix, autoArtRejParams, flagFig)

% patientCode = 'AMC062'; runCode = 1; suffix = 'HFA';

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

if nargin < 5
    flagFig = 1;
end
if nargin < 4
    autoArtRejParams = [7 .05]; % threshZ, threshShared
end

filename = sprintf('%s_preprocessed_ECoG/%s_TheWall1_run%i_preprocessed_%s.mat', workingDir, patientCode, runCode, suffix);
load(filename, 'ecog', 'patientInfo', 'stim32', 'CF32', 'stim128', 'CF128');

if strcmp(autoArtRejParams, 'manual')
    infoPatients = PF_infoPatients(patientCode);
    artifacts = infoPatients.artifacts{runCode};
    
    tb = patientInfo.time{1};
    [L, nElecs] = size(ecog);
    nArtifacts = length(artifacts);
    artifactsMat = zeros(L, nElecs);
    for idxArt = 1:nArtifacts
        artifactsMatTMP = zeros(L, nElecs);
        [~, artifactsSamp] = arrayfun(@(x) min(abs(tb-x)), artifacts{idxArt}(1:2));
        idxArtifacts = zeros(L, 1);
        idxArtifacts(artifactsSamp(1):artifactsSamp(2)-1) = 1;
        idxArtChan = ismember(1:nElecs, artifacts{idxArt}(3:end));
        artifactsMatTMP(:, idxArtChan) = repmat(idxArtifacts, 1, sum(idxArtChan));
        artifactsMat = artifactsMat + artifactsMatTMP;
    end
    artifacts = artifactsMat>0;
    rejRate = 100*sum(artifacts)/L;
    
    autoArtRejParams = [-1 -1];
else
    [artifacts, rejRate] = PF_getArtifacts(ecog, autoArtRejParams(1), autoArtRejParams(2), 20); % rejects 200 ms on each side of any sample which Z value exceeds 5; also, if more than 10% elecs share same artifact, then rej for all elecs
end

if flagFig > 0
    tb = patientInfo.time{1};
    nElecs = size(ecog, 2);
    ecogTMP = ecog;
    ecogTMP(:, patientInfo.noisyElecs{runCode}) = NaN;
    figure('Position', [6 47 2095 1209], 'Name', filename); subplot(4,1,1:3); imagesc(tb, 1:nElecs, (ecogTMP+artifacts)', [0 prctile(ecogTMP(:), 99.9)]);
    title(sprintf('threshZ = %.2g - threshAll = %.2g - mean rej rate = %.3g%%', autoArtRejParams, mean(rejRate)));
    set(gca, 'CLim', [0 1]);
    subplot(414); imagesc(tb, 1:nElecs, artifacts');

    if flagFig == 2
        while true
            idxElec = input('Which electrode do you want to observe? ([Enter] to continue)\n');
            if isempty(idxElec)
                break
            end
            figure; plot(tb, ecog(:, idxElec)); hold on;
            plot(tb(artifacts(:,idxElec)==1), ecog(artifacts(:,idxElec)==1, idxElec), 'r.');
            xlim(tb([1 end]));
        end
    end
end

patientInfo.rejRate = rejRate;
patientInfo.rejRateMean = mean(rejRate);

save(filename, 'ecog', 'patientInfo', 'stim32', 'CF32', 'stim128', 'CF128', 'artifacts');