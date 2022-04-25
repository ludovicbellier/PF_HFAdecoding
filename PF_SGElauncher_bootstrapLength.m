function PF_SGElauncher_bootstrapLength

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

pathSGE = [workingDir '_SGE/'];

nBootstrap = 100;
lengthTotal = 19072;
nElecsTotal = 347;
lengths = [15 30 60 90 120 150 180]; % in seconds

nSets = length(lengths);

nFbins = 32;
fBins = 1:nFbins;
for idxSet = 1:nSets
    lengthTMP = lengths(idxSet) * 100;
    memUsed = ceil((ceil(nElecsTotal*0.03 + 0.35) + 3) * lengthTMP / lengthTotal);
    begIdx = randperm(lengthTotal - lengthTMP) - 1;
    begIdx = begIdx(1:nBootstrap);
    for idxBS = 1:nBootstrap
        idxTotal = idxBS + (idxSet-1)*nBootstrap;
        
        fid = fopen(sprintf('%sparams%i.txt', pathSGE, idxTotal), 'w+');
        fprintf(fid, 'suffixPar = ''_bootstrapLength%i-%i''\n', idxSet, idxBS);
        fprintf(fid, 'nResamples = 5\n');
        fprintf(fid, 'flagTuneParam = 0\n');
        fprintf(fid, 'decodingLR = .001\n');
        fprintf(fid, 'idxSlice = [%i, %i]\n', [0 lengthTMP] + begIdx(idxBS));
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