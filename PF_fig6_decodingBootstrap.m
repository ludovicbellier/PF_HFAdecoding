function PF_fig6_decodingBootstrap

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

nBootstrap = 100;
nElecs = [5 10 20 40 80 160 320];
nSets = length(nElecs);

nFbins = 32;

fnameResults = sprintf('%sanalysis/supergrid_TheWall1_run1_HFA_bootstrap_results.mat', workingDir);
if exist(fnameResults, 'file') == 0
    superSaveR = zeros(nSets, nBootstrap, nFbins);
    superSaveR2 = zeros(nSets, nBootstrap, nFbins);
    superSaveTime = zeros(nSets, nBootstrap, nFbins);
    for idxSet = 1:nSets
        for idxBS = 1:nBootstrap
            for idxFreq = 1:nFbins
                try
                    load(sprintf('%sdecoding/supergrid_TheWall1_run1_HFA_bootstrap%i-%i_f%i.mat', workingDir, idxSet+10, idxBS, idxFreq), 'saveR', 'saveR2', 'totalTime');
                    superSaveR(idxSet, idxBS, idxFreq) = mean(saveR);
                    superSaveR2(idxSet, idxBS, idxFreq) = mean(saveR2);
                    superSaveTime(idxSet, idxBS, idxFreq) = totalTime;
                catch
                end
            end
        end
    end
    load(sprintf('%sdecoding/supergrid_TheWall1_run1_HFA_bootstrap11-1_f1.mat', workingDir), 'params');
    save(sprintf('%sanalysis/supergrid_TheWall1_run1_HFA_bootstrap_results.mat', workingDir), 'superSaveR', 'superSaveR2', 'superSaveTime', 'params', 'nElecs', 'nBootstrap', 'nFbins');
else
    load(fnameResults, 'superSaveR', 'nElecs', 'nBootstrap', 'nFbins');
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


%% Fig. 6A
figure('Position', [340 323 1277 837], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
values = squeeze(mean(values(1:nSets, :, :), 3, 'omitnan'));
mu = squeeze(mean(values, 2, 'omitnan'));
mu(isnan(mu)) = 0;
sem = squeeze(std(values, [], 2, 'omitnan')) ./ sqrt(nBootstrap);
errorbar(nElecs, mu, sem, 'Color', 'k', 'LineWidth', 1, 'CapSize', 3, 'Marker', '.');
xlabel('number of electrodes');
ylabel(yLabel);

mdl = fit(nElecs', mu, 'power2');
x = 1:347;
y = mdl(x);
hold on;
plot(x, y, 'Color', 'r', 'LineWidth', 1)
ylim([50 101]); box off; xlim([0 347]);

yThr = 80;
xThr = find(y>yThr, 1, 'first');
line([0 xThr], [yThr yThr], 'Color', 'r', 'LineStyle', ':');
line([xThr xThr], [0 yThr], 'Color', 'r', 'LineStyle', ':');
grid;
