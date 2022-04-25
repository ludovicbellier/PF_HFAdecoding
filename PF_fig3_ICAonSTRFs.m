global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end

%% define parameters
nFreq = 32; nLags = 75;
Ntotal = 2668;
inputArgList = {'AMC006', 1; 'AMC007', 1; 'AMC008', 1; 'AMC009', 1; 'AMC010', 1;...
    'AMC011', 1; 'AMC012', 2; 'AMC013', 1; 'AMC014', 1; 'AMC015', 1;...
    'AMC017', 1; 'AMC018', 2; 'AMC019', 1; 'AMC022', 1; 'AMC026', 1;...
    'AMC027', 1; 'AMC028', 1; 'AMC029', 3; 'AMC031', 1; 'AMC032', 1;...
    'AMC033', 1; 'AMC037', 1; 'AMC038', 1; 'AMC039', 2; 'AMC040', 2;...
    'AMC041', 2; 'AMC044', 1; 'AMC045', 1; 'AMC062', 1};
nCompInit = 10;

%% load STRF from all 29 patients
load(sprintf('%sanalysis/STRFmat_HFA_29pat.mat', workingDir));

%% center data
STRFmat = STRFmat - mean(STRFmat, 2);

%% perform ICA on the electrodes (to get several independent populations tuned to different spectrotemporal features)
Mdl = rica(STRFmat, nCompInit, 'VerbosityLevel', 1, 'IterationLimit', 5000); % STRFmat must be 2400x347
score_ICA = transform(Mdl, STRFmat);
coeff_ICA = Mdl.TransformWeights;

%% get explained variance for each component and sort components and show scree plot
explained = zeros(nCompInit, 1);
for idx = 1:nCompInit
    backProj = score_ICA(:,idx)*coeff_ICA(:,idx)';
    explained(idx) = 100-100*mean(var(STRFmat-backProj))/mean(var(STRFmat));
end
explained = sort(explained, 'descend');

thresh = 5;
figure; plot(explained, '.-'); hold on; line([0 nCompInit+1], [1 1]*thresh, 'LineStyle', '--', 'Color', 'r'); xlim([1 nCompInit]);
nComp = sum(explained >= thresh);

%% perform ICA on the electrodes (to get several independent populations tuned to different spectrotemporal features)
Mdl = rica(STRFmat, nComp, 'VerbosityLevel', 1, 'IterationLimit', 5000); % STRFmat must be 2400x347
score_ICA = transform(Mdl, STRFmat);
coeff_ICA = Mdl.TransformWeights;

% get explained variance for each component and sort components
explained = zeros(nComp, 1);
for idx = 1:nComp
    backProj = score_ICA(:,idx)*coeff_ICA(:,idx)';
    explained(idx) = 100-100*mean(var(STRFmat-backProj))/mean(var(STRFmat));
end
[explained, idxSort] = sort(explained, 'descend');
coeff_ICA = coeff_ICA(:, idxSort);
score_ICA = score_ICA(:, idxSort);

% flip components to have highest value as positive
for idx = 1:nComp
    if max(coeff_ICA(:, idx)) < abs(min(coeff_ICA(:, idx)))
        coeff_ICA(:, idx) = -1 .* coeff_ICA(:, idx);
        score_ICA(:, idx) = -1 .* score_ICA(:, idx);
    end
end


%%  Fig. 3B - Plot ICA components
nCompTMP = min([nComp 22]);
offset = 0;
comp_ICA_2D = zeros(nFreq, nLags, nComp);
figure('Position', [67 864 1367 356], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);
for idx = 1:nCompTMP
    idx2 = idx+offset;
    comp_ICA_2D(:, :, idx2) = reshape(score_ICA(:, idx2), nFreq, nLags);
    subplot(1, nCompTMP, idx); imagesc(comp_ICA_2D(:, :, idx2)); axis xy;
    set(gca, 'CLim', [-1 1].*max(abs(get(gca, 'CLim'))));
    title(sprintf('#%i - %.3g%% expl.', idx2, explained(idx2)));
end
colormap jet


%% Save ICA coefficients
save(sprintf('%sanalysis/STRFcoeffICA.mat', workingDir), 'coeff_ICA');


%% Save ICA components
save(sprintf('%sanalysis/ICAcomps.mat', workingDir), 'comp_ICA_2D');


%% Fig. 3C - Plot ICA component coefficients on MNI template brain
load(sprintf('%sanalysis/STRFmetrics_HFA_29pat.mat', workingDir));
idxSig = find(metrics(:, 7) == 1 & metrics(:, 6) > 0 & sum(metrics(:, 17:18), 2) == 0);

cMapRes = 512;
cMap = redblue(cMapRes);
cMap = [[0 0 0]; cMap(round(cMapRes/2)+1:end, :)];
for idx = 1:nComp
    values = nan(Ntotal, 1);
    values(idxSig) = coeff_ICA(:, idx).^.5; % compress values for visualization purpose
    values(values < 0) = 0;
    PF_plotValOnAnat(inputArgList, 'HFA', values, 'MNI', 1, 0, 'lh', [0 nan], 0, cMap);
    PF_plotValOnAnat(inputArgList, 'HFA', values, 'MNI', 1, 0, 'rh', [0 nan], 0, cMap);
end
