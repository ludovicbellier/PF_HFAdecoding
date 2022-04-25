function [allVal, idxOrig, metrics, metricsCols, elecCoord] = PF_plotValOnAnat(patientList, suffix, values, space, trigMask, figOutTrig, lateralFilter, valRange, centered, cmap)
%%% patientList = {'AMC026', 'AMC038', 'AMC062'}; space = 'MNI'; values = 'r'; trigMask = 2; highCutOff = []; centered = 0; strDir = '';
%%% choices for values => {numeric, '', 'sigR', 'sigR2', 'r', 'R2', 'iter', 'time', 'max', 'nan', 'rej'};
%%% patientList = 'AMC009'; suffix = 'ShortStories_70-150Hz';

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end
ftDir = '/home/knight/lbellier/DataWorkspace/_tools/git/fieldtrip/';
if exist('ft_defaults.m', 'file') == 0
    addpath(ftDir); ft_defaults;
end

% Define parameters
if nargin < 10
    cmap = 'default'; % 'default', 'jet', 'redblue', 'parula', 'lines', 'viridis'
end
if nargin < 9
    centered = 0;
end
if nargin < 8
    valRange = [nan nan];
end
if nargin < 7
    lateralFilter = '';
end
if nargin < 6
    figOutTrig = 0;
end
if nargin < 5
    trigMask = 0;
end
if nargin < 4
    space = 'patient';
end
if nargin < 3
    values = '';
end
elecSize = 12; % 30 % 50
edgeSize = 3; %2; % 2.5 % 4
labelId = 'fs';
toggleBadElecs = [1 0 0];%[1 0 0]; % display ref, noisy and epi elecs (0 -> no, 1 -> displayed, 2 -> displayed and tagged)
colorGrain = 0.01;
faceAlpha = 1;
flagCbar = 0;
flagVarSize = 0;
trigLegend = 0;
threshR2 = 0; % set to nan to deactivate this feature
anatDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/_anatomy/';


% Convert patientList to cell array when needed
if ~iscell(patientList)
    patientList = {patientList};
end
nPat = size(patientList, 1);


% Load electrode information
if exist(sprintf('%sanalysis/', workingDir), 'dir') == 0
    system(sprintf('mkdir %sanalysis', workingDir));
end
if exist(sprintf('%sanalysis/STRFs/', workingDir), 'dir') == 0
    system(sprintf('mkdir %sanalysis/STRFs', workingDir));
end
fname = sprintf('%sanalysis/STRFmetrics_%s_%ipat.mat', workingDir, suffix, nPat);
latDict = {'lh', 'rh', 'b'};
if nPat == 1 || ~exist(fname, 'file')
    [metrics, metricsCols] = PF_extractSTRFmetrics(patientList);
else
    load(fname, 'metrics', 'metricsCols');
end
nElecs = size(metrics, 1);
idxColSigR = ismember(metricsCols, 'sigR');
idxColR2 = ismember(metricsCols, 'R2');


% Load anatomical files
space = lower(space);
pathAnatTemplate = [anatDir 'AMC062/'];
if strcmp(space, 'patient')
    if nPat > 1
        space = 'mni';
        disp('Can''t use invidivual anatomy for several patients. MNI used instead.');
    elseif ~exist(sprintf('%s%s/freesurfer/', anatDir, patientList{1}), 'dir')
        space = 'mni';
        disp('Patient anatomy not available. MNI used instead.');
    end
end
switch space
    case 'patient'
        pathAnatTemplate = sprintf('%s%s/', anatDir, patientList{1});
        pial_lh = ft_read_headshape([pathAnatTemplate 'freesurfer/surf/lh.pial']);
        pial_rh = ft_read_headshape([pathAnatTemplate 'freesurfer/surf/rh.pial']);
        load([pathAnatTemplate patientList{1} '_elec_acpc_fr.mat'], 'elec_acpc_fr');
        elecCoord = ft_sortElec(elec_acpc_fr);
        suffixSpace = '';
    case 'fs'
        pial_lh = ft_read_headshape([pathAnatTemplate 'fsaverage/surf/lh.pial']);
        pial_rh = ft_read_headshape([pathAnatTemplate 'fsaverage/surf/rh.pial']);
        allElecs = cell(nPat, 1);
        for idxPat = 1:nPat
            pathAnatTemplate = sprintf('%s%s/', anatDir, patientList{idxPat});
            load([pathAnatTemplate patientList{idxPat} '_elec_fsavg_frs.mat'], 'elec_fsavg_frs');
            allElecs{idxPat} = ft_sortElec(elec_fsavg_frs);
        end
        elecCoord = ft_mergeElecSubsets(allElecs);
        suffixSpace = '-FS';
    case 'mni'
        pial = load([ftDir 'template/anatomy/surface_pial_left.mat']);
        pial_lh = pial.mesh;
        pial = load([ftDir 'template/anatomy/surface_pial_right.mat']);
        pial_rh = pial.mesh;
        clear pial mesh_lh mesh_rh
        allElecs = cell(nPat, 1);
        for idxPat = 1:nPat
            pathAnatTemplate = sprintf('%s%s/', anatDir, patientList{idxPat});
            load([pathAnatTemplate patientList{idxPat} '_elec_mni_frvr.mat'], 'elec_mni_frvr');
            allElecs{idxPat} = ft_sortElec(elec_mni_frvr);
        end
        elecCoord = ft_mergeElecSubsets(allElecs);
        suffixSpace = '-MNI';
end
elecCoord.label = sprintfc('%d', 1:nElecs);
idxUnk = find(~ismember(metrics(:, 3), [1 2]));
if ~isempty(idxUnk) % for bilateral cases, elecLat determined using x
    metrics(idxUnk, 3) = (elecCoord.elecpos(idxUnk, 1) > 0) + 1;
    save(fname, 'metrics', 'metricsCols');
end
allLat = unique(metrics(:, 3));
if length(allLat) == 1
    latImplant = latDict{allLat(1)};
else
    latImplant = lateralFilter;
end


% Filters
idxOrig = 1:nElecs;
idxFlagPlot = 1:nElecs;
labels = [];

if ~strcmp(lateralFilter, '')
    idxFlagPlot = idxFlagPlot(metrics(:, 3) == find(ismember(latDict, lateralFilter)));
    [idxFlagPlot, idxOrig, metrics, elecCoord, nElecs, values] = updateElecSet(idxFlagPlot, idxOrig, metrics, elecCoord, values, labels);
end

[refElec, noisyElecs, epilepticElecs] = getElecTags(metrics);
if isnumeric(values)
    idxSig = find(values ~= 0);
else
    if ~isnan(threshR2)
        idxSig = setdiff(find(metrics(idxFlagPlot, idxColSigR) == 1 & metrics(idxFlagPlot, idxColR2) > threshR2), [noisyElecs; epilepticElecs]);
    else
        idxSig = setdiff(find(metrics(idxFlagPlot, idxColSigR) == 1), [noisyElecs; epilepticElecs]);
    end
end
if trigMask == 2
    [idxFlagPlot, idxOrig, metrics, elecCoord, nElecs, values] = updateElecSet(idxSig, idxOrig, metrics, elecCoord, values, labels);    
    [refElec, noisyElecs, epilepticElecs] = getElecTags(metrics);
end


% Get values to plot
if isnumeric(values)
    strTitle = 'custom values';
    flagCbar = 1;
    [elecColors, colorMap, clim] = val2colorScale(values, colorGrain, valRange, centered, cmap);
    idxFlagPlot = setdiff(idxFlagPlot, find(isnan(values)));
    [idxFlagPlot, idxOrig, metrics, elecCoord, nElecs, values, labels, elecColors] = updateElecSet(idxFlagPlot, idxOrig, metrics, elecCoord, values, labels, elecColors);
    [refElec, noisyElecs, epilepticElecs] = getElecTags(metrics);
    allVal = values;
else
    strTitle = values;
    idxColVal = strcmp(metricsCols, values);
    allVal = metrics(:, idxColVal);
    switch values
        case ''
            elecColors = repmat([0 0 0], nElecs, 1);
            strTitle = 'coverage';
            
        case {'idPat', 'idElec'}
            flagCbar = 1;
            [elecColors, colorMap, clim] = val2colorScale(allVal, colorGrain, valRange, centered, cmap);
            
        case 'labels'
            trigLegend = 1;
            load(sprintf('%slabelColors.mat', workingDir), 'labelInfo_FS', 'labelInfo_MNI');
            labels = cell(nPat, 1);
            if strcmp(labelId, 'fs')
                for idxPat = 1:nPat
                    fnameLabels = sprintf('%s%s/%s_FsLabels.mat', anatDir, patientList{idxPat}, patientList{idxPat});
                    if exist(fnameLabels, 'file')
                        load(fnameLabels, 'FsLabels');
                        labels{idxPat} = FsLabels;
                    else
                        load(sprintf('%s%s/%s_labelTable.mat', workingDir, patientList{idxPat}, patientList{idxPat}), 'labelTable_FS', 'labelNames_FS');
                        [~, idxLabels] = max(labelTable_FS(:, 2:end), [], 2);
                        labels{idxPat} = labelNames_FS(idxLabels);
                    end
                end
                labelInfo = labelInfo_FS;
            elseif strcmp(labelId, 'mni')
                for idxPat = 1:nPat
                    load(sprintf('%s%s/%s_labelTable.mat', anatDir, patientList{idxPat}, patientList{idxPat}), 'labelTable_MNI', 'labelNames_MNI');
                    [~, idxLabels] = max(labelTable_MNI(:, 2:end), [], 2);
                    if length(idxLabels) ~= length(find(metrics(:,1)==idxPat))
                        idxLabels = idxLabels(1:length(find(metrics(:,1)==idxPat)));
                    end
                    labels{idxPat} = labelNames_MNI(idxLabels);
                end
                labelInfo = labelInfo_MNI;
            end
            labels = vertcat(labels{:});
            labels = labels(idxOrig);
            uniqueLabels = unique(labels);
            nLabels = length(uniqueLabels);
            idxColors = cellfun(@(x) find(strcmp(labelInfo(:, 1), x)), labels);
            idxUniqueColors = cellfun(@(x) find(strcmp(labelInfo(:, 1), x)), uniqueLabels);
            elecColors = cat(1, labelInfo{idxColors, 2});
            uniqueColors = cat(1, labelInfo{idxUniqueColors, 2});
            
        case {'sigR', 'sigR2', 'nanFail'} % binary values
            if strcmp(values, 'sigR') && trigMask ~= 2
                allVal = zeros(size(allVal));
                allVal(idxSig) = 1;
            end
            elecColors = repmat([1 .2 .2], nElecs, 1).*allVal;
            elecColors(allVal<0, :) = repmat([1 1 0], length(find(allVal<0)), 1);
            strTitle = sprintf('[%s] - %s - %i %s elec(s)', suffix, values, length(idxSig), values);
            
        case 'LR' % discrete values
            LRlist = unique(allVal);
            [~, allVal] = ismember(allVal, LRlist);
            flagCbar = 1;
            [elecColors, colorMap, clim] = val2colorScale(allVal, colorGrain, valRange, centered, 'hsv');

        otherwise % {'r', 'R2', 'iter', 'time', 'nanCount', 'maxCount'}
            if length(unique(allVal)) > 1
                [elecColors, colorMap, clim] = val2colorScale(allVal, colorGrain, valRange, centered, cmap);
                flagCbar = 1;
            else
                elecColors = repmat([0 0 0], nElecs, 1);
            end
            strTitle = sprintf('[%s] - %s - %i sig elec(s)', suffix, values, length(idxSig));
            if trigMask == 1
                idxNonSig = setdiff(1:nElecs, idxSig);
                elecColors(idxNonSig, :) = repmat([0 0 0], length(idxNonSig), 1);
            end
            values = allVal;
    end
end

if any(~isnan(valRange))
    idxFlagPlot = setdiff(idxFlagPlot, find(values<valRange(1) | values>valRange(2)));
	[idxFlagPlot, idxOrig, metrics, elecCoord, nElecs, values, labels, elecColors] = updateElecSet(idxFlagPlot, idxOrig, metrics, elecCoord, values, labels, elecColors);
    [refElec, noisyElecs, epilepticElecs] = getElecTags(metrics);
    strTitle = sprintf('%s - %i elecs displayed', strTitle, length(idxFlagPlot));
end

if any(toggleBadElecs < 1)
    badElecs = cell(3, 1);
    if toggleBadElecs(1) == 0 % ref
        badElecs{1} = refElec;
    end
    if toggleBadElecs(2) == 0 % noisy
        badElecs{2} = noisyElecs;
    end
    if toggleBadElecs(3) == 0 % epi
        badElecs{3} = epilepticElecs;
    end
    badElecs = cat(1, badElecs{:});
    if ~isempty(badElecs)
        idxFlagPlot = setdiff(idxFlagPlot, badElecs);
        [idxFlagPlot, idxOrig, metrics, elecCoord, nElecs, values, labels, elecColors] = updateElecSet(idxFlagPlot, idxOrig, metrics, elecCoord, values, labels, elecColors);
        [refElec, noisyElecs, epilepticElecs] = getElecTags(metrics);
    end
end

if flagVarSize == 1
    elecSize = round(metrics(:, 5) * 100) + 5;
end


% Plot values on surface
hf = figure('resize', 'off', 'Position', [46 276 1277 965], 'Color', [1 1 1], 'DefaultAxesFontSize', 5);%[1 1 640 640], 'Color', [1 1 1]);
if ~isempty(lateralFilter) && strcmpi('lh', lateralFilter)
    ft_plot_mesh(pial_lh, 'facealpha', faceAlpha);
    viewPoint = [-90 0];
elseif ~isempty(lateralFilter) && strcmpi('rh', lateralFilter)
    ft_plot_mesh(pial_rh, 'facealpha', faceAlpha);
    viewPoint = [90 0];
else
    ft_plot_mesh(pial_lh, 'facealpha', faceAlpha);
    ft_plot_mesh(pial_rh, 'facealpha', faceAlpha);
    if strcmp('lh', latImplant)
        viewPoint = [-90 0];
    else
        viewPoint = [90 0];
    end
end
if ~isempty(idxFlagPlot)
    try
        ft_plot_sens(elecCoord, 'elecsize', elecSize, 'facecolor', elecColors);
        fprintf('N_elecCoord=%i - N_elecColors=%i\n', length(elecCoord.label), length(elecColors));
    catch
        error('%i elecs available but trying to plot %i values', nElecs, size(elecColors, 1));
    end
    hold on;
    if toggleBadElecs(1) == 2
        h2 = scatter3(elecCoord.elecpos(refElec,1), elecCoord.elecpos(refElec,2), elecCoord.elecpos(refElec,3), elecSize(1)*edgeSize, [0 0.8 0], 'o', 'LineWidth', 2);
    end
    if toggleBadElecs(2) == 2
        h3 = scatter3(elecCoord.elecpos(noisyElecs,1), elecCoord.elecpos(noisyElecs,2), elecCoord.elecpos(noisyElecs,3), elecSize(1)*edgeSize, [0 0 0.8], 'o', 'LineWidth', 2);
    end
    if toggleBadElecs(3) == 2
        h4 = scatter3(elecCoord.elecpos(epilepticElecs,1), elecCoord.elecpos(epilepticElecs,2), elecCoord.elecpos(epilepticElecs,3), elecSize(1)*edgeSize, [0 0.6 0.8], 'o', 'LineWidth', 2);
    end
end
view([-90 0]); material dull; lighting gouraud; camlight;
view([90 0]); material dull; lighting gouraud; camlight;
view(viewPoint);

%%
% h = camlight('headlight');
% step = 0.25;
% speed = 0.1;
% for i = 0:step:36
%     if i>0
%         camorbit(10*step,0);
%         camlight(h,'headlight');
%     end
% %     pause(.01);
%     frame = getframe(hf);%, [110 80 370 280]);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if i == 0
%         imwrite(imind,cm,'allElecsOnSurface_FS.gif','gif','DelayTime',speed*step,'Loopcount',inf);
%     else
%         imwrite(imind,cm,'allElecsOnSurface_FS.gif','gif','DelayTime',speed*step,'WriteMode','append');
%     end
% end
%%

elecLabels = cell(nElecs, 1);
for idxElec = 1:nElecs
    elecLabels{idxElec} = sprintf('%s e%i', patientList{metrics(idxElec, 1)}, metrics(idxElec, 2));
end
dcm = datacursormode(hf);
set(dcm, 'update', {@PF_cursorOutputFunction, elecLabels, [], values}, 'Interpreter', 'none');

patientListStr = sprintf('%s, ', patientList{:, 1});
patientListStr = ['[' patientListStr(1:end-2) ']'];
if nPat <= 5
    title(sprintf('%s - %s', patientListStr, strTitle), 'Interpreter', 'none');
else
    title(sprintf('%i patients - %s', nPat, strTitle), 'Interpreter', 'none');
end
if flagCbar == 1
    hCbar = colorbar('Location', 'south'); colormap(colorMap); set(gca, 'CLim', double(clim));
    if exist('LRlist', 'var')
        hCbar.TickLabels = LRlist;
    end
end
if trigLegend == 1
    h5 = plot(nan(2, nLabels), '.', 'MarkerSize', 20);
    set(h5, {'color'}, mat2cell(uniqueColors, ones(nLabels, 1)));
    hLegend = [h2; h3; h4; h5];
    labelLegend = ['ref elec'; 'noisy elecs'; 'epileptic elecs'; uniqueLabels];
    toggleLegend = false(3, 1);
    for idxLegend = 1:3
        if isempty(eval(sprintf('h%i.XData', idxLegend+1)))
            toggleLegend(idxLegend) = true;
        end
    end
    hLegend(toggleLegend) = [];
    labelLegend(toggleLegend) = [];
    legend(hLegend, labelLegend, 'Interpreter', 'none', 'Location', 'EastOutside');
    set(dcm, 'update', {@PF_cursorOutputFunction, elecLabels, labels}, 'Interpreter', 'none');
end
if figOutTrig == 1
    set(gcf, 'Color', [1 1 1]);
    pause(1);
    if nPat == 1
        saveas(gcf, sprintf('%sanalysis/STRFs/%s_%s_valOnAnat%s.png', workingDir, patientList{1}, suffix, suffixSpace));
    else
        saveas(gcf, sprintf('%sanalysis/STRFs/%ipatients_%s_STRFvalsOnAnat%s-%s.png', workingDir, nPat, suffix, suffixSpace, latImplant));
    end
end

    function [idxFlagPlot, idxOrig, metrics, elecCoord, nElecs, values, labels, elecColors] = updateElecSet(idxFlagPlot, idxOrig, metrics, elecCoord, values, labels, elecColors)
        elecCoord = ft_selectElecSubset(elecCoord, idxFlagPlot);
        metrics = metrics(idxFlagPlot, :);
        if isnumeric(values)
            values = values(idxFlagPlot);
        end
        if exist('elecColors', 'var')
            elecColors = elecColors(idxFlagPlot, :);
        end
        if ~isempty(labels)
            labels = labels(idxFlagPlot);
        end
        nElecs = length(idxFlagPlot);
        idxOrig = idxOrig(idxFlagPlot);
        idxFlagPlot = 1:nElecs;
    end

    function [refElec, noisyElecs, epilepticElecs] = getElecTags(metrics)
        refElec = find(metrics(:, end-2));
        noisyElecs = find(metrics(:, end-1));
        epilepticElecs = find(metrics(:, end));
    end
end