function PF_SGElauncher_modeling(patientCode, runCode, modelType, suffix, elecs, paramFile)

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end
pathSGE = [workingDir '_SGE/'];

if ~exist('paramFile', 'var')
    paramFile = [];
end
if ~exist('elecs', 'var') || isempty(elecs)
    try
        patientInfo = PF_infoPatients(patientCode);
        assert(~isempty(patientInfo));
    catch
        load(sprintf('%s_preprocessed_ECoG/%s_TheWall1_run%i_preprocessed_%s.mat', workingDir, patientCode, runCode, suffix), 'patientInfo');
    end
    elecs = patientInfo.ecogElec{1};
end
nElecs = length(elecs);
memUsed = 3;

if any(strcmp(modelType, {'decoding', 'recon'}))
    sigElecs = elecs;
    nSigElecs = length(sigElecs);
    memUsed = ceil(nSigElecs*0.03 + 0.35) + 3;
    switch modelType
        case 'decoding'
            nElecs = 32;
%         case 'recon'
%             nElecs = 128;
    end
%     elecs = 1:nElecs;
end

cellOut = cell(nElecs, 1);
for i = 1:nElecs
    cellOut{i} = sprintf('PF-%s-%i-%i-%s-%i-%s', patientCode, runCode, elecs(i), modelType, paramFile, suffix);
end

fid = fopen([pathSGE, 'SGE_joblist.txt'], 'w+');
for i = 1:length(cellOut)
    fprintf(fid, '%s\n', cellOut{i});
end
fclose(fid);

% memUsed = 15;
fid = fopen([pathSGE, 'SGE_options.txt'], 'w+');
% fprintf(fid, '-M ludovic.bellier@berkeley.edu -m eas -l mem_free=%iG', memUsed);
fprintf(fid, '-l mem_free=%iG', memUsed);
fclose(fid);

pause(2);
system(sprintf('submit -s %sPF_SGE_runLinRegESE.sh -f %sSGE_joblist.txt -o %sSGE_options.txt', pathSGE, pathSGE, pathSGE));