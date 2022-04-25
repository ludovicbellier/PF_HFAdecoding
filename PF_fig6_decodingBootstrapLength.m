function PF_fig6_decodingBootstrapLength

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

nBootstrap = 100;
lengths = [15 30 60 90 120 150 180]; % in seconds
nSets = length(lengths);

nFbins = 32;

fnameResults = sprintf('%sanalysis/supergrid_TheWall1_run1_HFA_bootstrapLength_results.mat', workingDir);
if exist(fnameResults, 'file') == 0
    superSaveR = zeros(nSets, nBootstrap, nFbins);
    superSaveR2 = zeros(nSets, nBootstrap, nFbins);
    superSaveTime = zeros(nSets, nBootstrap, nFbins);
    for idxSet = 1:nSets
        for idxBS = 1:nBootstrap
            for idxFreq = 1:nFbins
                try
                    load(sprintf('%sdecoding/supergrid_TheWall1_run1_HFA_bootstrapLength%i-%i_f%i.mat', workingDir, idxSet, idxBS, idxFreq), 'saveR', 'saveR2', 'totalTime');
                    superSaveR(idxSet, idxBS, idxFreq) = mean(saveR);
                    superSaveR2(idxSet, idxBS, idxFreq) = mean(saveR2);
                    superSaveTime(idxSet, idxBS, idxFreq) = totalTime;
                catch
                end
            end
        end
    end
    load(sprintf('%sdecoding/supergrid_TheWall1_run1_HFA_bootstrapLength1-1_f1.mat', workingDir), 'params');
    save(sprintf('%sanalysis/supergrid_TheWall1_run1_HFA_bootstrapLength_results.mat', workingDir), 'superSaveR', 'superSaveR2', 'superSaveTime', 'params', 'lengths', 'nBootstrap', 'nFbins');
else
    load(fnameResults, 'superSaveR', 'lengths', 'nBootstrap', 'nFbins');
end


values = superSaveR;
values(values==0) = nan;

% normalize based on best prediction accuracy using all 347 sig elecs
bestR = zeros(100, nFbins);
for idxFreq = 1:nFbins
    load(sprintf('%sdecoding/supergrid_TheWall1_run1_HFA_f%i.mat', workingDir, idxFreq), 'save_r');
    bestR(:, idxFreq) = save_r;
end
bestRmu = mean(bestR, 'omitnan');

values = bsxfun(@rdivide, values, permute(bestRmu, [1 3 2])) * 100;
yLabel = '% of best Pearson''r';


%% Fig. 6B
figure('Position', [340 323 1277 837], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
values = squeeze(mean(values(1:nSets, :, :), 3, 'omitnan'));
mu = squeeze(mean(values, 2, 'omitnan'));
mu(isnan(mu)) = 0;
sem = squeeze(std(values, [], 2, 'omitnan')) ./ sqrt(nBootstrap);
errorbar(lengths, mu, sem, 'Color', 'k', 'LineWidth', 1, 'CapSize', 3, 'Marker', '.');
xlabel('dataset duration');
ylabel(yLabel);

mdl = fit(lengths', mu, 'power2');
x = 1:191;
y = mdl(x);
hold on;
plot(x, y, 'Color', 'r', 'LineWidth', 1)
ylim([50 101]); box off; xlim([0 191]);

yThr = 90;
xThr = find(y>yThr, 1, 'first');
line([0 xThr], [yThr yThr], 'Color', 'r', 'LineStyle', ':', 'LineWidth', .5);
line([xThr xThr], [0 yThr], 'Color', 'r', 'LineStyle', ':', 'LineWidth', .5);
grid;
