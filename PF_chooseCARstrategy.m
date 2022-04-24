function elecSubsetsOut = PF_chooseCARstrategy(data, noisyElecs, goodElecs, freqBand, fname, elecSubsets)

threshGaussGOF = .95;
batchSize = 16;
flagNaive = 0;
rerefMethod = 'median';
nElecs = size(data.trial{1}, 1);

% covariance matrix for raw data
cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfreq = [[60 120 180]-.5; [60 120 180]+.5]';
cfg.bpfilter = 'yes';
cfg.bpfreq = freqBand;
dataTMP = ft_preprocessing(cfg, data);
covRaw = getCovarianceMatrix(data.trial{1}', noisyElecs);
[xRaw, yRaw, aR2Raw, yRaw2, aR2Raw2, covRawZ] = fitCovHistGaussian(covRaw);

% covariance matrix for CAR on all electrodes
cfg = [];
cfg.reref = 'yes';
cfg.refchannel = goodElecs;
cfg.refmethod = rerefMethod;
dataCARall = ft_preprocessing(cfg, dataTMP);
covCARall = getCovarianceMatrix(dataCARall.trial{1}', noisyElecs);
[xCARall, yCARall, aR2CARall, yCARall2, aR2CARall2, covCARallZ] = fitCovHistGaussian(covCARall);

% covariance matrix for CAR on electrode subsets
if ~exist('elecSubsets', 'var')
    nElecBatches = nElecs/batchSize;
    if nElecBatches - floor(nElecBatches) >= .5
        nElecBatches = ceil(nElecBatches);
    else
        nElecBatches = floor(nElecBatches);
    end
    batchIndices = cell(nElecBatches, 1);
    for idxBatch = 1:nElecBatches
        if idxBatch < nElecBatches
            batchIndices{idxBatch} = 1+(idxBatch-1)*batchSize:idxBatch*batchSize;
        else
            batchIndices{idxBatch} = 1+(idxBatch-1)*batchSize:nElecs;
        end
    end
    batchID = zeros(nElecBatches);
    for idxBatch = 1:nElecBatches
        rowMeans = mean(covCARall(batchIndices{idxBatch}, :), 'omitnan');
        batchMeans = cellfun(@(x) mean(rowMeans(x), 'omitnan'), batchIndices);
        batchIdAll = zeros(100, nElecBatches);
        for idxIter = 1:100
            batchIdTMP= kmeans(batchMeans, 2);
            if batchIdTMP(1) == 1 && mean(batchMeans(batchIdTMP==1)) > mean(batchMeans(batchIdTMP==2))
                batchIdTMP = 2./batchIdTMP;
            end
            batchIdAll(idxIter, :) = batchIdTMP;
        end
        batchID(idxBatch, :) = mode(batchIdAll);
    end
    batchID = round(mean(batchID));
    if batchID(1) == 2
        batchID = 2./batchID;
    end
    if length(unique(batchID)) > 1
        elecSubsets = {[batchIndices{batchID==1}], [batchIndices{batchID==2}]};
    else
        flagNaive = 1;
        idx1 = 2^nextpow2(nElecs/2);
        idx2 = 2^(nextpow2(nElecs/2)-1);
        if abs(2*idx2 - nElecs) > abs(2*idx1 - nElecs)
            idxEndSub1 = idx1;
        else
            idxEndSub1 = idx2;
        end
        elecSubsets = {1:idxEndSub1, idxEndSub1+1:nElecs};
    end
end
cfg = [];
cfg.reref = 'yes';
cfg.refmethod = rerefMethod;
dataSubsets = zeros(size(dataCARall.trial{1}));
for idxSubset = 1:length(elecSubsets)
    cfg.refchannel = elecSubsets{idxSubset}(ismember(elecSubsets{idxSubset}, goodElecs));
    dataPartialCAR = ft_preprocessing(cfg, dataTMP);
    dataSubsets(elecSubsets{idxSubset}, :) = dataPartialCAR.trial{1}(elecSubsets{idxSubset}, :);
end
dataCARsub = dataCARall;
dataCARsub.trial{1} = dataSubsets;
covCARsub = getCovarianceMatrix(dataCARsub.trial{1}', noisyElecs);
[xCARsub, yCARsub, aR2CARsub, yCARsub2, aR2CARsub2, covCARsubZ] = fitCovHistGaussian(covCARsub);

% select best strategy
if aR2CARall < threshGaussGOF && aR2CARsub > aR2CARall && flagNaive == 0
    CARallColor = 'k';
    CARsubColor = 'r';
    elecSubsetsOut = elecSubsets;
else
    CARallColor = 'r';
    if flagNaive == 1
        CARsubColor = 'b';
    else
        CARsubColor = 'k';
    end
    elecSubsetsOut = {};
end

% plot covariance matrices, histograms and fitted gaussians
figure('Position', [46 613 1200 620]);
subplot(231); imagesc(covRaw, getRobustClim(covRaw)); colorbar; title('Raw');
subplot(232); imagesc(covCARall, getRobustClim(covCARall)); colorbar; title(sprintf('CARall (%s)', rerefMethod), 'Color', CARallColor);
subplot(233); imagesc(covCARsub, getRobustClim(covCARsub)); colorbar; title(sprintf('CARsub\n%s', getSubsetsString(elecSubsets)), 'Color', CARsubColor);
subplot(234); histogram(covRawZ(:), 50); hold on; plot(xRaw, yRaw, 'r', 'LineWidth', 2);
plot(xRaw, yRaw2, 'g', 'LineWidth', 2);
title(sprintf('gauss_aR2=%.3f - gauss2_aR2=%.3f', aR2Raw, aR2Raw2), 'Interpreter', 'none', 'Color', 'k');
subplot(235); histogram(covCARallZ(:), 50); hold on; plot(xCARall, yCARall, 'r', 'LineWidth', 2);
plot(xCARall, yCARall2, 'g', 'LineWidth', 2);
title(sprintf('gauss_aR2=%.3f - gauss2_aR2=%.3f', aR2CARall, aR2CARall2), 'Interpreter', 'none', 'Color', CARallColor);
subplot(236); histogram(covCARsubZ(:), 50); hold on; plot(xCARsub, yCARsub, 'r', 'LineWidth', 2);
plot(xCARsub, yCARsub2, 'g', 'LineWidth', 2);
title(sprintf('gauss_aR2=%.3f - gauss2_aR2=%.3f', aR2CARsub, aR2CARsub2), 'Interpreter', 'none', 'Color', CARsubColor);

% save figure
pause(2);
saveas(gcf, fname);
% pause(2); close all;

% declare nested functions
    function M_out = getCovarianceMatrix(M_in, noisyElecs)
        M_out = cov(M_in);
        M_out(noisyElecs, :) = nan; M_out(:, noisyElecs) = nan;
        M_out(logical(eye(size(M_out)))) = nan;
    end

    function [x, yFit_Gauss1, aR2_Gauss1, yFit_Gauss2, aR2_Gauss2, covMat] = fitCovHistGaussian(covMat)
        okFlag = false;
        margin = .5;
        
        while ~okFlag
            try
                covMat(covMat<=prctile(covMat(:), margin)) = NaN;
                covMat(covMat>=prctile(covMat(:), 100-margin)) = NaN;
                
                [y, x] = histcounts(covMat(:), 50);
                x = x(1:50);
                [fitresult, gof] = fit(x', y', 'gauss1');
                aR2_Gauss1 = gof.adjrsquare;
                yFit_Gauss1 = fitresult(x);
                
                [fitresult, gof] = fit(x', y', 'gauss2');
                aR2_Gauss2 = gof.adjrsquare;
                yFit_Gauss2 = fitresult(x);
                
                okFlag = true;
            catch
                margin = margin + .5;
            end
        end
    end

    function robClim = getRobustClim(M_in, prctiles)
        if nargin < 2
            prctiles = [10 90];
        end
        robClim = prctile(M_in(:), prctiles);
    end

    function [strSubsets] = getSubsetsString(chanIdx)
        strSubsets = sprintf('{');
        L = size(chanIdx, 2);
        for i = 1:L
            if length(unique(diff(chanIdx{i}))) == 1
                strSubsets = [strSubsets, sprintf('%i:%i', chanIdx{i}([1 end]))]; %#ok
            else
                strSubsets = [strSubsets, '[']; %#ok
                jumpIdx = find(diff(chanIdx{i}) > 1);
                jumpIdxMat = [chanIdx{i}([1 jumpIdx+1]); chanIdx{i}([jumpIdx end])]';
                L2 = size(jumpIdxMat, 1);
                for j = 1:L2
                    strSubsets = [strSubsets, sprintf('%i:%i', jumpIdxMat(j, [1 2]))]; %#ok
                    if j < L2
                        strSubsets = [strSubsets, ', ']; %#ok
                    end
                end
                strSubsets = [strSubsets, ']']; %#ok
            end
            if i < L
                strSubsets = [strSubsets, sprintf(', ')]; %#ok
            end
        end
        strSubsets = [strSubsets, sprintf('}')];
    end
end