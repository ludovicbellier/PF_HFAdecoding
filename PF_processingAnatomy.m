function PF_processingAnatomy(patientCode)
% patientCode = 'AMC007';

global workingDir;
if isempty(workingDir)
    workingDir = '/home/knight/lbellier/DataWorkspace/_projects/PinkFloyd/';
end
ftDir = '/home/knight/lbellier/DataWorkspace/_tools/git/fieldtrip/';
if exist('ft_defaults.m', 'file') == 0
    addpath(ftDir); ft_defaults;
end

%% Anatomical workflow step 1
patientInfo = PF_infoPatients(patientCode);
anatDir = [workingDir '_anatomy/' patientCode '/'];
if exist(anatDir, 'dir') == 0
    system(sprintf('mkdir %s', anatDir));
end
fshome = '/usr/local/freesurfer_x86_64-6.0.0/';
ft_hastoolbox('spm12', 1);

figPos = [1 42 1530 646];
% figPos = [1 -560 1530 1248];
% figPos = [1998 143 1260 1100];

if ~isempty(patientInfo.MR)
    %% Step 2
    mri = ft_read_mri(patientInfo.MR);
    
    %% Step 3
    ft_determine_coordsys(mri);
    close all
    % X is the left-right axis, and + is on the right => left-to-right
    % orientation
    
    %% Step 4
    cfg           = [];
    cfg.method    = 'interactive';
    cfg.coordsys  = 'acpc';
    mri_acpc = ft_volumerealign(cfg, mri);
    
    %% Step 5
    cfg           = [];
    cfg.filename  = [anatDir patientCode '_MR_acpc'];
    cfg.filetype  = 'nifti';
    cfg.parameter = 'anatomy';
    ft_volumewrite(cfg, mri_acpc);
    
    %% Step 6
    subdir     = anatDir;
    mrfile     = [anatDir patientCode '_MR_acpc.nii'];
    system(['export FREESURFER_HOME=' fshome '; ' ...
        'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
        'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...
        'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all'])
    
    %% Step 7
    pial_lh = ft_read_headshape([anatDir 'freesurfer/surf/lh.pial']);
    pial_lh.coordsys = 'acpc';
    figure;
    subplot(1,2,1);
    ft_plot_mesh(pial_lh);
    view([270 0]); % left view
    material dull;
    lighting gouraud;
    camlight;
    
    pial_rh = ft_read_headshape([anatDir 'freesurfer/surf/rh.pial']);
    pial_rh.coordsys = 'acpc';
    subplot(1,2,2);
    ft_plot_mesh(pial_rh);
    view([90 0]); % right view
    material dull;
    lighting gouraud;
    camlight;
    
    %% Step 8
    fsmri_acpc = ft_read_mri([anatDir 'freesurfer/mri/T1.mgz']);
    fsmri_acpc.coordsys = 'acpc';
else
    fsmri_acpc = ft_read_mri([ftDir 'template/anatomy/single_subj_T1_1mm.nii']);
    fsmri_acpc.coordsys = 'acpc';
end

%% Step 9
ct = ft_read_mri(patientInfo.CT);

%% Step 10
ft_determine_coordsys(ct);
close all

%% Step 11
cfg           = [];
cfg.method    = 'interactive';
cfg.coordsys  = 'ctf';
ct_ctf = ft_volumerealign(cfg, ct);

%% Step 12
ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc');
% load([pathAnat patientCode '_CT_acpc']); % AMC028, once FS done

%% Step 13
cfg             = [];
cfg.method      = 'spm';
cfg.spmversion  = 'spm12';
cfg.coordsys    = 'acpc';
cfg.viewresult  = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);

%% Step 15
cfg           = [];
cfg.filename  = [anatDir patientCode '_CT_acpc_f'];
cfg.filetype  = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);

%% Step 16
hdr = ft_read_header(patientInfo.filename{1});
hdr = hdr.label; %#ok
save([anatDir patientCode '_hdrLabels.mat'], 'hdr'); % for use on local computer

%% Step 17 - performed on a local computer
% cfg         = [];
% cfg.channel = hdr.label;
% elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);

%% Step 19
if ~isempty(patientInfo.MR)
    fsmri_acpc = ft_read_mri([anatDir 'freesurfer/mri/T1.mgz']);
    fsmri_acpc.coordsys = 'acpc';
else
    fsmri_acpc = ft_read_mri([ftDir 'template/anatomy/single_subj_T1_1mm.nii']);
    fsmri_acpc.coordsys = 'acpc';
end
load([anatDir patientCode '_elecCoord.mat'])
ft_plot_ortho(fsmri_acpc.anatomy, 'transform', fsmri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');

%% Step 20
save([anatDir patientCode '_elec_acpc_f.mat'], 'elec_acpc_f');
delete ([anatDir patientCode '_elecCoord.mat'], [anatDir patientCode '_hdrLabels.mat']);
if ~isempty(patientInfo.MR)
    %% Step 21
    cfg           = [];
    cfg.method    = 'cortexhull';
    cfg.fshome    = fshome;
    cfg.headshape = [anatDir 'freesurfer/surf/lh.pial'];
%     cfg.smooth_steps = 60;
    cfg.outer_surface_sphere = 80;
    hull_lh = ft_prepare_mesh(cfg);
%     hull_lh = inflateHull(hull_lh, ft_read_headshape([pathAnat 'freesurfer/surf/lh.pial']), []);
    cfg.headshape = [anatDir 'freesurfer/surf/rh.pial'];
    hull_rh = ft_prepare_mesh(cfg);
%     hull_rh = inflateHull(hull_rh, ft_read_headshape([pathAnat 'freesurfer/surf/rh.pial']), []);

    %% Step 22
    save([anatDir patientCode '_hull_lh.mat'], 'hull_lh');
    save([anatDir patientCode '_hull_rh.mat'], 'hull_rh');
end

%% Step 23
if ~isempty(patientInfo.MR)
    load([anatDir patientCode '_hull_lh.mat']);
    load([anatDir patientCode '_hull_rh.mat']);
else
    load('/home/knight/ecog_tmp/Recon_Functions/template_hull.mat');
    clear mesh mesh_lh mesh_rh
end
latImplant = patientInfo.laterality;
if strcmp(latImplant, 'lh')
    hull = hull_lh;
    viewPoint = [270 0];
else
    hull = hull_rh;
    viewPoint = [90 0];
end

load([anatDir patientCode '_elec_acpc_f.mat']);
grids = patientInfo.ecogGrids;
strips = patientInfo.ecogStrips;
gridsNstrips = [grids strips];
gridsNstripsID = [ones(1, length(grids)) 2*ones(1, length(strips))];
warpMethods = {'dykstra2012', 'hermes2010'};

elec_acpc_fr = elec_acpc_f;
for g = 1:numel(gridsNstrips)
    tic
    [elecTMP, idxTMP] = selectElecSubset(elec_acpc_f, gridsNstrips{g});
    
    cfg             = [];
    cfg.keepchannel = 'yes';
    cfg.elec        = elecTMP;
    cfg.method      = 'headshape';
    cfg.headshape   = hull;
    cfg.warp        = warpMethods{gridsNstripsID(g)};
    cfg.feedback    = 'yes';
    elecTMP2 = ft_electroderealign(cfg);
%     deformArr = [0.01 0.1 0.5 1 2];
%     clear elecTMP2
%     for i = 1:length(deformArr)
%         cfg.deformweight = deformArr(i);
%         elecTMP2{i} = ft_electroderealign(cfg);
%     end
%     cfg.warp = 'hermes2010';
%     elecTMP2{i+1} = ft_electroderealign(cfg);
%     for i = 1:length(deformArr)+1
%         figure('Position', figPos);
%         ft_plot_mesh(hull, 'facealpha', 0.8); ft_plot_sens(elecTMP, 'elecsize', 15, 'facecolor', 'k'); ft_plot_sens(elecTMP2{i}, 'elecsize', 15, 'facecolor', 'r');
%         if i < length(deformArr)+1
%             title(sprintf('dykstra2012 (deforweight=%.3g) - normDist=%.3g', deformArr(i), norm(elecTMP2{i}.elecpos - elecTMP.elecpos)/size(elecTMP.elecpos, 1)));
%         else
%             title(sprintf('hermes2010 - normDist=%.3g', norm(elecTMP2{i}.elecpos - elecTMP.elecpos)/size(elecTMP.elecpos, 1)));
%         end
%         view(viewPoint); material dull; lighting gouraud; camlight;
%     end
%     elecTMP2 = elecTMP2{4};
    
%     figure('Position', figPos);
%     ft_plot_mesh(hull, 'facealpha', 0.8); ft_plot_sens(elecTMP, 'elecsize', 15, 'facecolor', 'k'); ft_plot_sens(elecTMP2, 'elecsize', 15, 'facecolor', 'r');
%     title(sprintf('%s - normDist=%.3g', warpMethods{gridsNstripsID(g)}, norm(elecTMP2.elecpos - elecTMP.elecpos)/size(elecTMP.elecpos, 1)));
%     view(viewPoint); material dull; lighting gouraud; camlight;
    
    elec_acpc_fr.chanpos(idxTMP, :) = elecTMP2.chanpos;
    elec_acpc_fr.elecpos = elec_acpc_fr.chanpos;
    toc
end

figure('Position', figPos);
ft_plot_mesh(hull, 'facealpha', 0.8); ft_plot_sens(elec_acpc_f, 'elecsize', 15, 'facecolor', 'k'); ft_plot_sens(elec_acpc_fr, 'elecsize', 15, 'facecolor', 'r');
title(sprintf('normDist=%.3g', norm(elec_acpc_fr.elecpos - elec_acpc_f.elecpos)/size(elec_acpc_f.elecpos, 1)));
view(viewPoint); material dull; lighting gouraud; camlight;

%% Step 24
% figure;
% ft_plot_mesh(pial); ft_plot_mesh(hull, 'facecolor', 'r');
% ft_plot_sens(elec_acpc_fr);
% view(viewPoint); material dull; lighting gouraud; camlight;
% [hull2, elec_acpc_fr2] = inflateHull(hull, pial, elec_acpc_fr);
% figure;
% ft_plot_mesh(pial); ft_plot_mesh(hull2, 'facecolor', 'r');
% ft_plot_sens(elec_acpc_fr2);
% view(viewPoint); material dull; lighting gouraud; camlight;

% correct bad position
% load([pathAnat patientCode '_elec_acpc_fr.mat'], 'elec_acpc_fr');
% hull = hull_lh_r80;
% pial_lh = ft_read_headshape([pathAnat 'freesurfer/surf/lh.pial']);
% pial = pial_lh;
% viewPoint = [-90 0];
% figure;
% % ft_plot_mesh(hull); ft_plot_sens(elec_acpc_fr_her, 'label', 'label');
% ft_plot_mesh(pial); ft_plot_mesh(hull, 'FaceAlpha', .4); ft_plot_sens(elec_acpc_fr, 'label', 'label');
% view(viewPoint); material dull; lighting gouraud; camlight;
% elecTMP2 = elec_acpc_fr;
% elecLabel = '22';
% idxLabel = strcmp(elecTMP2.label, elecLabel);
% elecTMP2.elecpos(idxLabel, :) = cursor_info.Position;
% elecTMP2.chanpos(idxLabel, :) = cursor_info.Position;

%% Step 25
save([anatDir patientCode '_elec_acpc_fr.mat'], 'elec_acpc_fr');

%% Step 26 - Volume-based normalization -> MNI template
if ~isempty(patientInfo.MR)
    fsmri_acpc = ft_read_mri([anatDir 'freesurfer/mri/T1.mgz']);
    fsmri_acpc.coordsys = 'acpc';
    cfg            = [];
    cfg.nonlinear  = 'yes';
    cfg.spmversion = 'spm12';
    fsmri_mni = ft_volumenormalise(cfg, fsmri_acpc);
end

%% Step 27
elec_mni_frv = elec_acpc_fr;
elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni.params, elec_acpc_fr.elecpos, 'individual2sn');
elec_mni_frv.chanpos = ft_warp_apply(fsmri_mni.params, elec_acpc_fr.chanpos, 'individual2sn');
elec_mni_frv.coordsys = 'mni';
save([anatDir patientCode '_elec_mni_frv.mat'], 'elec_mni_frv');

%% Step 27.5
load('/home/knight/ecog/ecog_scripts/Recon_Functions/template_hull.mat');
% now in /home/knight/ecog/ecog_scripts/data_analysis/Recon_Functions/template_hull.mat
clear mesh mesh_lh mesh_rh
latImplant = patientInfo.laterality;
if strcmp(latImplant, 'lh')
    hull = hull_lh;
    viewPoint = [270 0];
else
    hull = hull_rh;
    viewPoint = [90 0];
end

load([anatDir patientCode '_elec_mni_frv.mat']);
grids = patientInfo.ecogGrids;
strips = patientInfo.ecogStrips;
gridsNstrips = [grids strips];
gridsNstripsID = [ones(1, length(grids)) 2*ones(1, length(strips))]; % 2*ones(1, length(gridsNstrips));
warpMethods = {'dykstra2012', 'hermes2010'};

elec_mni_frvr = elec_mni_frv;
for g = 1:numel(gridsNstrips)
    tic
    [elecTMP, idxTMP] = ft_selectElecSubset(elec_mni_frv, gridsNstrips{g});
    
    cfg             = [];
    cfg.keepchannel = 'yes';
    cfg.elec        = elecTMP;
    cfg.method      = 'headshape';
    cfg.headshape   = hull;
    cfg.warp        = warpMethods{gridsNstripsID(g)};
    cfg.feedback    = 'yes';
    elecTMP2 = ft_electroderealign(cfg);
    deformArr = [0.01 0.1 0.5 1 2];
    clear elecTMP2
    for i = 1:length(deformArr)
        cfg.deformweight = deformArr(i);
        elecTMP2{i} = ft_electroderealign(cfg);
    end
    cfg.warp = 'hermes2010';
    elecTMP2{i+1} = ft_electroderealign(cfg);
    for i = 1:length(deformArr)+1
        figure('Position', figPos);
        ft_plot_mesh(hull, 'facealpha', 0.8); ft_plot_sens(elecTMP, 'elecsize', 15, 'facecolor', 'k'); ft_plot_sens(elecTMP2{i}, 'elecsize', 15, 'facecolor', 'r');
        if i < length(deformArr)+1
            title(sprintf('dykstra2012 (deforweight=%.3g) - normDist=%.3g', deformArr(i), norm(elecTMP2{i}.elecpos - elecTMP.elecpos)/size(elecTMP.elecpos, 1)));
        else
            title(sprintf('hermes2010 - normDist=%.3g', norm(elecTMP2{i}.elecpos - elecTMP.elecpos)/size(elecTMP.elecpos, 1)));
        end
        view(viewPoint); material dull; lighting gouraud; camlight;
    end
    %     elecTMP2 = elecTMP2{6};
    
    elec_mni_frvr.chanpos(idxTMP, :) = elecTMP2.chanpos;
    elec_mni_frvr.elecpos = elec_mni_frvr.chanpos;
    toc
end

figure('Position', figPos);
ft_plot_mesh(hull, 'facealpha', 0.8); ft_plot_sens(elec_mni_frv, 'elecsize', 15, 'facecolor', 'k'); ft_plot_sens(elec_mni_frvr, 'elecsize', 15, 'facecolor', 'r');
title(sprintf('normDist=%.3g', norm(elec_mni_frvr.elecpos - elec_mni_frv.elecpos)/size(elec_mni_frv.elecpos, 1)));
view(viewPoint); material dull; lighting gouraud; camlight;

% %% Step 28
% elecsize = 8;
% figure;
% subplot(2,2,1); ft_plot_mesh(hull, 'facealpha', 0.9); ft_plot_sens(elec_mni_frv, 'elecsize', elecsize);
% view(viewPoint); material dull; lighting gouraud; camlight;
% subplot(2,2,2); ft_plot_mesh(hull, 'facealpha', 0.9); ft_plot_sens(elec_mni_frvr, 'elecsize', elecsize);
% view(viewPoint); material dull; lighting gouraud; camlight;
% subplot(2,2,3); ft_plot_mesh(pial, 'facealpha', 0.9); ft_plot_sens(elec_mni_frv, 'elecsize', elecsize);
% view(viewPoint); material dull; lighting gouraud; camlight;
% subplot(2,2,4); ft_plot_mesh(pial, 'facealpha', 0.9); ft_plot_sens(elec_mni_frvr, 'elecsize', elecsize);%, 'facecolor', 'r');
% view(viewPoint); material dull; lighting gouraud; camlight;

%% Step 29
save([anatDir patientCode '_elec_mni_frvr.mat'], 'elec_mni_frvr');

%% Step 30 - Surface-based normalization -> Freesurfer average template
latImplant = patientInfo.laterality;
cfg           = [];
cfg.channel   = 'all';
cfg.elec      = elec_acpc_fr;
cfg.method    = 'headshape';
cfg.headshape = sprintf('%sfreesurfer/surf/%s.pial', anatDir, latImplant);
cfg.warp      = 'fsaverage';
cfg.fshome    = fshome;
elec_fsavg_frs = ft_electroderealign(cfg);

pial_FS = ft_read_headshape(sprintf('%s/subjects/fsaverage/surf/%s.pial', fshome, latImplant));
pial_FS.coordsys = 'fsaverage';
figure('Position', figPos);
ft_plot_mesh(pial_FS, 'facealpha', 0.8); ft_plot_sens(elec_fsavg_frs, 'elecsize', 15, 'facecolor', 'k');
view(viewPoint); material dull; lighting gouraud; camlight;

% %% Step 31
% if strcmp(latImplant, 'lh')
%     viewPoint = [270 0];
% else
%     viewPoint = [90 0];
% end
% fspial = ft_read_headshape(sprintf('%s/subjects/fsaverage/surf/%s.pial', fshome, latImplant));
% fspial.coordsys = 'fsaverage';
% ft_plot_mesh(fspial);
% ft_plot_sens(elec_fsavg_frs);
% view(viewPoint); material dull; lighting gouraud; camlight;

%% Step 32
save([anatDir patientCode '_elec_fsavg_frs.mat'], 'elec_fsavg_frs');

%% Step 33
atlas_MNI = ft_read_atlas([ftDir 'template/atlas/aal/ROI_MNI_V4.nii']);
atlas_MNI.coordsys = 'mni';
atlas_FS = ft_read_atlas([anatDir 'freesurfer/mri/aparc+aseg.mgz']);
atlas_FS.coordsys = 'acpc';
mri_MNI = ft_read_mri([ftDir 'template/anatomy/single_subj_T1_1mm.nii']);
mri_MNI.coordsys = 'mni';

%% Step 34
load([anatDir patientCode '_elec_acpc_fr.mat']);
load([anatDir patientCode '_elec_mni_frvr.mat']);

cfg = [];
cfg.parameter = 'tissue';
atlas_MNIinterp = ft_sourceinterpolate(cfg, atlas_MNI, mri_MNI);
atlas_MNIinterp.coordsys = 'mni';

nElecs = length(elec_acpc_fr.label);
labelTable = cell(nElecs, 3);
labelTable(:, 1) = elec_acpc_fr.label;
nElecs = length(elec_mni_frvr.label);
labelTable = cell(nElecs, 118);
labelTable(:, 1) = elec_mni_frvr.label;

% tic
cfg = [];
cfg.roi = elec_mni_frvr.chanpos;
cfg.atlas = atlas_MNI;
cfg.inputcoord = 'mni';
cfg.output = 'multiple';
cfg.minqueryrange = 5;
cfg.maxqueryrange = 29;
labels_MNI = ft_volumelookup(cfg, atlas_MNIinterp);

cfg.roi = elec_acpc_fr.chanpos;
cfg.atlas = atlas_FS;
cfg.inputcoord = 'acpc';
labels_FS = ft_volumelookup(cfg, atlas_FS);

labCount = [labels_MNI.count]';
labelTable_MNI = [str2double(elec_mni_frvr.label), labCount];%round(labCount ./ sum(labCount, 2) * 100)];
labelNames_MNI = labels_MNI(1).name;

labCount = [labels_FS.count]';
labelTable_FS = [str2double(elec_acpc_fr.label), labCount];%round(labCount ./ sum(labCount, 2) * 100)];
labelNames_FS = labels_FS(1).name;

% for idxElec = 1:nElecs
% %     [~, idx] = max(labels_MNI(idxElec).count);
% %     labelTable(idxElec, 2) = labels_MNI(idxElec).name(idx);
%     idx = find(labels_MNI(idxElec).count>0);
%     if length(idx) == 1
%         labelTable(idxElec, 2) = labels_MNI(idxElec).name(idx);
%     else
%         percent = round(labels_MNI(idxElec).count(idx)/sum(labels_MNI(idxElec).count(idx))*100);
%         labNames = labels_MNI(idxElec).name(idx);
%         for idxLab = 1:length(labNames)
%             labNames{idxLab} = sprintf('%s (%i%%)', labNames{idxLab}, percent(idxLab));
%         end
%         [~, idx2] = sort(percent, 'descend');
%         labTMP = sprintf('%s - ', labNames{idx2});
%         labelTable{idxElec, 2} = labTMP(1:end-3);
%     end
% %     [~, idx] = max(labels_FS(idxElec).count);
% %     labelTable(idxElec, 3) = labels_FS(idxElec).name(idx);
% end
% toc

save([anatDir patientCode '_labelTable.mat'], 'labels_MNI', 'labelTable_MNI', 'labelNames_MNI', 'labels_FS', 'labelTable_FS', 'labelNames_FS');
% save([pathAnat patientCode '_labelTable.mat'], 'labels_MNI', 'labelTable_MNI', 'labelNames_MNI');
% save([pathAnat patientCode '_labelTable_FS.mat'], 'labels_FS', 'labelTable_FS', 'labelNames_FS');

% %% Step 35 - Plot everything
% elecsize = 15;
% % Patient
% if ~isempty(patientInfo.MR)
%     pial_lh = ft_read_headshape([pathAnat 'freesurfer/surf/lh.pial']);
%     pial_rh = ft_read_headshape([pathAnat 'freesurfer/surf/rh.pial']);
%     load([pathAnat patientCode '_hull_lh.mat']);
%     load([pathAnat patientCode '_hull_rh.mat']);
% else
%     pial = load([ftDir 'template/anatomy/surface_pial_left.mat']);
%     pial_lh = pial.mesh;
%     pial = load([ftDir 'template/anatomy/surface_pial_right.mat']);
%     pial_rh = pial.mesh;
%     load('/home/knight/ecog_tmp/Recon_Functions/template_hull.mat');
%     clear pial mesh_lh mesh_rh
% end
% latImplant = patientInfo.laterality;
% if strcmp(latImplant, 'lh')
%     hull_pat = hull_lh;
%     pial_pat = pial_lh;
%     viewPoint = [270 0];
% else
%     hull_pat = hull_rh;
%     pial_pat = pial_rh;
%     viewPoint = [90 0];
% end
% load([pathAnat patientCode '_elec_acpc_f.mat']);
% load([pathAnat patientCode '_elec_acpc_fr.mat']);
% 
% % figure;
% % ft_plot_mesh(pial_pat); ft_plot_mesh(hull_pat, 'facecolor', 'r');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % 
% % figure;
% % subplot(2,2,1); ft_plot_mesh(hull_pat, 'facealpha', 0.9); ft_plot_sens(elec_acpc_f, 'elecsize', elecsize);%, 'label', 'label');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % title('Without brain-shift compensation');
% % subplot(2,2,2); ft_plot_mesh(hull_pat, 'facealpha', 0.9); ft_plot_sens(elec_acpc_fr, 'elecsize', elecsize);%, 'label', 'label');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % title('With brain-shift compensation - Dykstra 2012');
% % subplot(2,2,3); ft_plot_mesh(pial_pat, 'facealpha', 0.9); ft_plot_sens(elec_acpc_f, 'elecsize', elecsize);%, 'label', 'label');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % subplot(2,2,4); ft_plot_mesh(pial_pat, 'facealpha', 0.9); ft_plot_sens(elec_acpc_fr, 'elecsize', elecsize);%, 'label', 'label');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % 
% % figure;
% % ft_plot_mesh(pial_pat); ft_plot_sens(elec_acpc_fr, 'elecsize', 15);%, 'label', 'label');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% 
% % MNI
% pial = load([ftDir 'template/anatomy/surface_pial_left.mat']);
% pial_lh = pial.mesh;
% pial = load([ftDir 'template/anatomy/surface_pial_right.mat']);
% pial_rh = pial.mesh;
% load('/home/knight/ecog_tmp/Recon_Functions/template_hull.mat');
% clear pial mesh_lh mesh_rh
% latImplant = patientInfo.laterality;
% if strcmp(latImplant, 'lh')
%     hull_MNI = hull_lh;
%     pial_MNI = pial_lh;
%     viewPoint = [270 0];
% else
%     hull_MNI = hull_rh;
%     pial_MNI = pial_rh;
%     viewPoint = [90 0];
% end
% load([pathAnat patientCode '_elec_mni_frv.mat']);
% load([pathAnat patientCode '_elec_mni_frvr.mat']);
% 
% % figure;
% % subplot(2,2,1); ft_plot_mesh(hull_MNI, 'facealpha', 0.9); ft_plot_sens(elec_mni_frv, 'elecsize', elecsize);
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % subplot(2,2,2); ft_plot_mesh(hull_MNI, 'facealpha', 0.9); ft_plot_sens(elec_mni_frvr, 'elecsize', elecsize);
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % subplot(2,2,3); ft_plot_mesh(pial_MNI, 'facealpha', 0.9); ft_plot_sens(elec_mni_frv, 'elecsize', elecsize);
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % subplot(2,2,4); ft_plot_mesh(pial_MNI, 'facealpha', 0.9); ft_plot_sens(elec_mni_frvr, 'elecsize', elecsize);
% % view(viewPoint); material dull; lighting gouraud; camlight;
% 
% % FSavg
% pial_FS = ft_read_headshape(sprintf('%s/subjects/fsaverage/surf/%s.pial', fshome, latImplant));
% pial_FS.coordsys = 'fsaverage';
% load([pathAnat patientCode '_elec_fsavg_frs.mat']);
% 
% % figure;
% % ft_plot_mesh(pial_FS);
% % ft_plot_sens(elec_fsavg_frs);
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % 
% % % all
% % figure;
% % subplot(1,3,1); ft_plot_mesh(pial_pat); ft_plot_sens(elec_acpc_fr, 'elecsize', elecsize); title('patient''s brain');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % subplot(1,3,2); ft_plot_mesh(pial_MNI); ft_plot_sens(elec_mni_frvr, 'elecsize', elecsize); title('MNI template');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% % subplot(1,3,3); ft_plot_mesh(pial_FS); ft_plot_sens(elec_fsavg_frs, 'elecsize', elecsize); title('Freesurfer template');
% % view(viewPoint); material dull; lighting gouraud; camlight;
% 
% % with anatomical labels
% % ColorOrder = get(gca, 'ColorOrder');
% % close(gcf);
% colOrdTMP = jet;
% idx = round(linspace(1, 64, 20));
% ColorOrder = colOrdTMP(idx, :);
% % ColorOrder(idx,:) = colOrdTMP;
% colOrdTMP1 = lines;
% colOrdTMP2 = prism;
% ColorOrder = [colOrdTMP1(1:7,:); colOrdTMP2(1:6,:); ColorOrder];
% 
% 
% load([pathAnat patientCode '_labelTableOrig.mat']);
% labels = labelTable(:, 3);
% uniqueLabels = unique(labels);
% figure('Position', figPos);
% ft_plot_mesh(pial_pat, 'facealpha', 0.8);
% for idxLabel = 1:length(uniqueLabels)
%     elecTMP = selectElecSubset(elec_acpc_fr, str2double(labelTable(strcmp(labels, uniqueLabels(idxLabel)), 1)));
%     ft_plot_sens(elecTMP, 'elecsize', elecsize, 'facecolor', ColorOrder(idxLabel, :));
% end
% view(viewPoint); material dull; lighting gouraud; camlight;
% legend(['patient_brain'; uniqueLabels], 'Interpreter', 'none');
% title('patient''s brain');
% 
% labels = labelTable(:, 3);
% uniqueLabels = unique(labels);
% figure('Position', figPos);
% ft_plot_mesh(pial_MNI, 'facealpha', 0.8);
% for idxLabel = 1:length(uniqueLabels)
%     elecTMP = selectElecSubset(elec_mni_frvr, str2double(labelTable(strcmp(labels, uniqueLabels(idxLabel)), 1)));
%     ft_plot_sens(elecTMP, 'elecsize', elecsize, 'facecolor', ColorOrder(idxLabel, :), 'label', 'label');
% end
% view(viewPoint); material dull; lighting gouraud; camlight;
% legend(['MNI_brain'; uniqueLabels], 'Interpreter', 'none');
% title('MNI brain');
% 
% labels = labelTable(:, 3);
% uniqueLabels = unique(labels);
% figure('Position', figPos);
% ft_plot_mesh(pial_FS, 'facealpha', 0.8);
% for idxLabel = 1:length(uniqueLabels)
%     elecTMP = selectElecSubset(elec_fsavg_frs, str2double(labelTable(strcmp(labels, uniqueLabels(idxLabel)), 1)));
%     ft_plot_sens(elecTMP, 'elecsize', elecsize, 'facecolor', ColorOrder(idxLabel, :));
% end
% view(viewPoint); material dull; lighting gouraud; camlight;
% legend(['MNI_brain'; uniqueLabels], 'Interpreter', 'none');
% title('FSavg brain');
% 
% %%
% % fname = 'ribbon';
% % fname = 'aparc+aseg';
% % fname = 'aparc.a2009s+aseg';
% % fnameIn = [pathAnat 'freesurfer/mri/' fname '.mgz'];
% % fnameOut = [pathAnat fname '.nii'];
% % fshome     = '/usr/local/freesurfer_x86_64-6.0.0/';
% % system(sprintf('export FREESURFER_HOME=%s; source $FREESURFER_HOME/SetUpFreeSurfer.sh; mri_convert --in_type mgz --out_type nii --out_orientation RAS %s %s', fshome, fnameIn, fnameOut));