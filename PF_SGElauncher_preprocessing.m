function PF_SGElauncher_preprocessing(patientCode, runCode)

pathSGE = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/_SGE/';

patientInfo = PF_infoPatients(patientCode);
nElecs = length(patientInfo.ecogElec{1});
memUsed = ceil(nElecs/16);

paramList = {'HFA'};
nPar = length(paramList);

fid = fopen([pathSGE, 'SGE_joblist.txt'], 'w+');
for idxPar = 1:nPar
    fprintf(fid, '%s-%i-%s\n', patientCode, runCode, paramList{idxPar});
end
fclose(fid);

fid = fopen([pathSGE, 'SGE_options.txt'], 'w+');
% fprintf(fid, '-M ludovic.bellier@berkeley.edu -m eas -l mem_free=%iG', memUsed);
fprintf(fid, '-l mem_free=%iG', memUsed);
fclose(fid);

pause(2);
system(sprintf('submit -s %sPF_SGE_runPreprocECoG.sh -f %sSGE_joblist.txt -o %sSGE_options.txt', pathSGE, pathSGE, pathSGE));