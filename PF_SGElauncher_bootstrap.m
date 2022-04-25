function PF_SGElauncher_bootstrap

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

pathSGE = [workingDir '_SGE/'];

nBootstrap = 100;
nElecsTotal = 347;
nElecs = [5 10 20 40 80 160 320];

nSets = length(nElecs);

nFbins = 32;
fBins = 1:nFbins;
for idxSet = 5:nSets
    nTMP = nElecs(idxSet);
    memUsed = ceil(nTMP*0.03 + 0.35) + 3;
    for idxBS = 1:nBootstrap
        idxTotal = idxBS + (idxSet-1)*nBootstrap;
        
        elecIdx = randperm(nElecsTotal)-1;
        elecIdx = sort(elecIdx(1:nTMP));
        elecIdxStr = strjoin(sprintfc('%i, ', elecIdx));
        
        fid = fopen(sprintf('%sparams%i.txt', pathSGE, idxTotal), 'w+');
        fprintf(fid, 'suffixPar = ''_bootstrap%i-%i''\n', idxSet+10, idxBS);
        fprintf(fid, 'nResamples = 5\n');
        fprintf(fid, 'flagTuneParam = 0\n');
        fprintf(fid, 'decodingLR = .001\n');
        fprintf(fid, 'elecIdx = [%s]\n', elecIdxStr(1:end-2));
        fclose(fid);
        
        pause(2);
        
        cellOut = cell(nFbins, 1);
        for i = 1:nFbins
            cellOut{i} = sprintf('PF-supergrid-1-%i-decoding-%i-HFA', fBins(i), idxTotal);
        end
        
        fid = fopen([pathSGE, 'SGE_joblist.txt'], 'w+');
        for i = 1:length(cellOut)
            fprintf(fid, '%s\n', cellOut{i});
        end
        fclose(fid);
        
        fid = fopen([pathSGE, 'SGE_options.txt'], 'w+');
        % fprintf(fid, '-M ludovic.bellier@berkeley.edu -m eas -l mem_free=%iG', memUsed);
        fprintf(fid, '-l mem_free=%iG', memUsed);
        fclose(fid);
        
        pause(2);
        system(sprintf('submit -s %sPF_SGE_runLinRegESE.sh -f %sSGE_joblist.txt -o %sSGE_options.txt', pathSGE, pathSGE, pathSGE));
    end
end