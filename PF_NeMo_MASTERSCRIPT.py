"""
Script to compute encoding and decoding models using NeMo functions.

In the Pink Floyd | HFA decoding study, we assess the relationship between the auditory spectrogram of a Pink Floyd song
 used as a stimulus and the elicited High-Frequency Activity from 29 patients with pharmaco-resistant epilepsy.

This script loads the preprocessed data, prepare it for modeling by computing the feature lag matrix, and computes
 either encoding or decoding models, with many controllable parameters.
"""

import sys
import time

import numpy as np
import scipy.io as sio

sys.path.extend(['/home/knight/lbellier/DataWorkspace/_tools/git/NeMo/'])
import nemo as nm


########################################################################################################################
#  Initialization
########################################################################################################################
# Get input parameters (when called through Sun Grid Engine)
# argin = 'PF-AMC062-1-195-encoding--HFB'
argin = sys.argv[1]

# Initialize parameters and variables
params = {'flag_fig': 0,
          'flag_tune': 1,
          'random_state': None,
          'fs': 100,
          'stim_type': 'stim32',
          'lags_encoding': (.75, 0),
          'lr_encoding': .001,
          'lags_decoding': (.75, 0),
          'lr_decoding': .001,
          'split_ratio': (.7, .2, .1),
          'group_length': 2,
          'fixed_test': (),
          'vocal_boundaries': (),
          'flag_log': 0,
          'flag_artifact_correct': 1,
          'slice_idx': (),
          'n_splits': 100,
          'hidden_layer_sizes': (64, 64,),
          'alpha': .01,
          'l1_ratio': .5,
          'n_best_of_n': 3,
          'criterion_best_of_n': 'mse_val',
          'batch_size': 200,
          'early_stopping': True,
          'n_iter_no_change': 10,
          'tol': 0,
          'max_iter': 200,
          'min_iter': 5,
          'max_splits': 200,
          'lim_nan': 50,
          'lim_max': 200,
          'algo': 'rMLRwES',
          'scaler_type': 'robust',
          'scaler_iqr': (2, 98),
          'n_cv_tuning': 5,
          'verbose': 1}

dataset_id, patient_code, run_code, target_idx, params['model_type'], param_file = argin.split('-')[:6]
run_code = int(run_code)
target_idx = int(target_idx) - 1
suffix = '-'.join(argin.split('-')[6:])
if suffix != '':
    suffix = '_' + suffix
elec_idx = []
target_sets = []
classes = []
cv_results = ()
hyperparam_grid = {}
best_hyperparam = []
working_dir = '/home/knight/lbellier/DataWorkspace/_projects/'
data_dir = '_preprocessed_ECoG/'

# Load data
time_start = time.time()
if dataset_id in ['PF', 'PinkFloyd']:
    working_dir += 'PinkFloyd/'
    stim_name = 'TheWall1_run%d' % run_code
    params['vocal_boundaries'] = (13, 80)  # time boundaries of the vocals (sec), for stratified split
    params['flag_log'] = 1
    if params['model_type'] == 'recon':
        params['fixed_test'] = (60, 75)  # fixed test set to be reconstructed (sec)
else:
    stim_name = ''
fname_in = '%s%s_%s_preprocessed%s.mat' % (data_dir, patient_code, stim_name, suffix)
data = sio.loadmat(working_dir + fname_in)
if params['model_type'] == 'recon':
    params['stim_type'] = 'stim128'
    params['algo'] = 'MLP'
if target_idx < 0:  # multi-output case
    target_idx = np.arange(data[params['stim_type']].shape[1])

# (optional) Load custom parameters from an external text file
if param_file != '':
    param_suffix = '_params' + param_file
    param_fname = working_dir + '_SGE/params' + param_file + '.txt'
    f = open(param_fname, 'r')
    param_update = [x for x in f.read().split('\n')[:-1]]
    f.close()
    param_update = [x for x in param_update if x]
    param_update = ['params[\'' + x[0:x.find(' =')] + '\'] ' + x[x.find('= '):] for x in param_update]
    [exec(x) for x in param_update]
else:
    param_suffix = ''


########################################################################################################################
#  Data preparation
########################################################################################################################
# Create classes for stratified split
if dataset_id in ['PF', 'PinkFloyd', 'CR', 'continuousRec', 'TIMIT', 'AL']:
    classes = nm.create_classes(data[params['stim_type']], **params)
elif dataset_id in ['SS', 'MVS', 'ShortStories']:
    classes = np.asarray(data['strata'][:, 0], dtype='float64') - 1

# Define features and target
if params['model_type'] == 'encoding':
    features = data[params['stim_type']]
    if params['flag_log'] == 1:
        if features.min() <= 0:
            features += sys.float_info.epsilon
        features = np.log(features)
    target = data['ecog'][:, target_idx]
    if 'artifacts' in data.keys():
        artifacts = data['artifacts'][:, target_idx]
    else:
        artifacts = np.zeros(target.shape, dtype='uint8')
    target_prefix = 'e'
    params['lags'] = params['lags_encoding'][0]
    params['offset'] = params['lags_encoding'][1]
    params['learning_rate'] = params['lr_encoding']
elif params['model_type'] in ['decoding', 'recon']:
    if len(elec_idx) == 0 and 'supergrid' not in patient_code:
        try:
            STRF_metrics = sio.loadmat(working_dir + 'analysis/STRFmetrics_HFB_29pat.mat')
            idx_pat = np.where(STRF_metrics['patList'][:, 0] == patient_code)[0][0] + 1
            STRF_metrics = STRF_metrics['metrics'][STRF_metrics['metrics'][:, 0] == idx_pat, :]
            elec_idx = np.intersect1d(np.intersect1d(np.where(STRF_metrics[:, 6] == 1)[0],
                                                     np.where(STRF_metrics[:, 5] >= .01)[0]),
                                      np.where(np.sum(STRF_metrics[:, -2:], 1) == 0)[0])
            if elec_idx[-1] > data['ecog'].shape[1]:
                elec_idx = np.arange(0, data['ecog'].shape[1])
        except FileNotFoundError:
            print('no STRFmetrics*.mat file - using all electrodes for decoding/recon instead')
            elec_idx = np.arange(0, data['ecog'].shape[1])
    else:
        elec_idx = np.arange(0, data['ecog'].shape[1])
    features = data['ecog'][:, elec_idx]
    target = data[params['stim_type']][:, target_idx]
    if params['flag_log'] == 1:
        if target.min() <= 0:
            target += .1  # sys.float_info.epsilon
        target = np.log(target)
    if 'artifacts' in data.keys() and params['flag_artifact_correct'] == 1:
        artifacts = data['artifacts'][:, elec_idx]
    else:
        artifacts = np.zeros(features.shape, dtype='uint8')
    [features, target, artifacts, classes] = nm.fix_artifacts(features, target, artifacts, classes,
                                                              'replace', **params)
    target_prefix = 'f'
    params['lags'] = params['lags_decoding'][0]
    params['offset'] = params['lags_decoding'][1]
    params['learning_rate'] = params['lr_decoding']
else:
    sys.exit("Wrong model type specified: please choose between 'encoding', 'decoding' and 'recon'.")

# (optional) Take a subset of the data
# e.g., params['slice_idx'] = [10000, 19071]
if len(params['slice_idx']) > 0:
    slice_idx = slice(params['slice_idx'][0], params['slice_idx'][1])
    features = features[slice_idx, :]
    target = target[slice_idx]
    artifacts = artifacts[slice_idx]
    classes = classes[slice_idx]

# Build the feature lag matrix
[features, target, artifacts, classes] = nm.build_lag_matrix(features, target, artifacts, classes, **params)
if params['model_type'] == 'encoding':
    [features, target, artifacts, classes] = nm.fix_artifacts(features, target, artifacts, classes,
                                                              'remove', **params)
params['n_feat'] = features.shape[1]
params['n_feat_y'] = int(params['n_feat'] / (params['lags'] * params['fs']))
params['n_target'] = target.shape[1]

# (optional) Estimate feature auto-correlation span (use this to find a proper value for chunk duration)
# nm.get_feat_corr(features, 10, 2, params['fs'])

# Create groups for split
groups = nm.create_groups(classes, **params)

# (optional) Plot preprocessed (modeling-ready) data
# nm.plot_prepared_data(features, target, classes, groups, **params)


########################################################################################################################
#  Hyperparameter tuning
########################################################################################################################
if params['flag_tune'] > 0:
    if params['algo'] == 'rMLRwES':
        hyperparam_grid = {'learning_rate': np.logspace(-3, 2, 11)}
    elif params['algo'] == 'ridge':
        hyperparam_grid = {'alpha': np.logspace(-2, 2, 9)}
    elif params['algo'] == 'lasso':
        hyperparam_grid = {'alpha': np.logspace(-5, 1, 13)}
    elif params['algo'] == 'elasticNet':
        hyperparam_grid = {'alpha': np.logspace(-5, 1, 13), 'l1_ratio': [.1, .25, .5, .75, .9]}
    elif params['algo'] == 'MLP':
        hyperparam_grid = {'alpha': np.logspace(-2, 2, 9)}
    else:
        hyperparam_grid = {}

    [best_hyperparam, cv_results, y_pred_cv] = nm.grid_search(features, target, groups, classes, params,
                                                              hyperparam_grid)
    params.update(best_hyperparam)

    if params['flag_tune'] == 2:
        fname_out_cv = '%s/%s_%s%s%s_%s%d-tune%s_results.mat' % (params['model_type'], patient_code, stim_name, suffix,
                                                                 param_suffix, target_prefix, target_idx + 1,
                                                                 params['algo'])
        if params['random_state'] is None:
            params['random_state'] = 'None'
        if params['flag_log'] == 1:
            y_pred_cv = np.exp(y_pred_cv)
        mdict = {'cv_results': cv_results, 'hyperparam_grid': hyperparam_grid, 'y_pred_cv': y_pred_cv, 'params': params}
        sio.savemat(working_dir + fname_out_cv, mdict=mdict)
        sys.exit()


########################################################################################################################
#  Model fitting
########################################################################################################################
# Initialize saving variables and loop over resamples
if isinstance(target_idx, int):
    fname_out = '%s/%s_%s%s%s_%s%d.mat' % (params['model_type'], patient_code, stim_name, suffix, param_suffix,
                                           target_prefix, target_idx + 1)
else:
    fname_out = '%s/%s_%s%s%s_all%s.mat' % (params['model_type'], patient_code, stim_name, suffix, param_suffix,
                                            target_prefix.upper())
save_coefs = np.zeros((params['n_splits'], params['n_feat_y'], int(params['lags'] * params['fs'])))
save_iter = np.zeros((params['n_splits']))
save_runtime = np.zeros((params['n_splits']))
save_metrics = []
save_split_indices = []
save_split_ratios = []
save_y_pred = []
idx = 0
max_iter_counter = 0
nan_counter = 0
flag_last_chance = 0
while idx < params['n_splits']:
    print('Processing split #{}/{}'.format(idx + 1, params['n_splits']))
    split_indices = nm.stratified_group_shuffle_split(groups, classes, **params)
    [feature_sets, target_sets] = nm.perform_split(features, target, split_indices, **params)
    [feature_sets, scaler] = nm.scale(feature_sets, **params)
    if params['algo'] == 'rMLRwES':
        mdl = nm.RobustMultipleLinearRegressionEarlyStopping(**params)
        if params['early_stopping']:
            mdl.fit(feature_sets[0], target_sets[0], feature_sets[1], target_sets[1])
        else:
            mdl.fit(np.concatenate(feature_sets[:2]), np.concatenate(target_sets[:2]))

    elif params['algo'] in ['ridge', 'lasso', 'elasticNet']:
        mdl = nm.RegularizedRegressionCustomEstimator(**params)
        mdl.fit(np.concatenate(feature_sets[:2]), np.concatenate(target_sets[:2]))

    elif params['algo'] == 'MLP':
        params['criterion_best_of_n'] = 'r2_test'
        mdl = nm.MLPRegressorCustomEarlyStopping(**params)
        mdl, _ = nm.best_of_n(mdl, feature_sets, target_sets, scaler, params)

    else:
        mdl = []

    [metrics, save_coefs[idx, :, :], yPred, save_iter[idx], save_runtime[idx], _] = \
        nm.get_model_output(mdl, feature_sets, target_sets, scaler, **params)
    if mdl.is_fitted_ and not np.isnan(metrics[2, 0]):
        save_metrics.append(metrics)
        save_split_indices.append(split_indices)
        save_split_ratios.append([len(split_indices[i]) / features.shape[0] for i in range(len(split_indices))])
        save_y_pred.append(yPred)
        idx += 1
    elif save_iter[idx] == params['max_splits']:
        max_iter_counter += 1
        print('Warning: max iteration cap has been reached')
    else:
        nan_counter += 1
        print("Warning: nan score - model didn't converge")
    if nan_counter == params['lim_nan']:
        if flag_last_chance == 0:
            if params['flag_tune'] > 0:
                idx_new_hp = np.where(hyperparam_grid['learning_rate'] == best_hyperparam['learning_rate'])[0][0] - 1
                if idx_new_hp >= 0:
                    new_learning_rate = hyperparam_grid['learning_rate'][idx_new_hp]
                else:
                    break
            else:
                new_learning_rate = params['learning_rate'] / 10
            params['learning_rate'] = new_learning_rate
            save_metrics = []
            save_split_indices = []
            save_split_ratios = []
            save_y_pred = []
            idx = 0
            nan_counter = 0
            flag_last_chance = 1
        else:
            break
    if max_iter_counter == params['lim_max']:
        break
save_r = np.asarray([x[2, 0] for x in save_metrics])
save_r2 = np.asarray([x[2, 1] for x in save_metrics])
save_mse = np.asarray([x[2, 2] for x in save_metrics])


########################################################################################################################
#  Model visualization and saving
########################################################################################################################
# Plot output
if params['flag_fig'] > 1:
    if params['model_type'] in ['encoding', 'decoding']:
        nm.plot_strf(save_coefs, save_r, save_r2, **params)
    if params['model_type'] in ['decoding', 'recon'] and len(params['fixed_test']) > 0:
        nm.plot_y_pred(save_y_pred, target_sets[-1], **params)

# Save modeling results
if params['random_state'] is None:
    params['random_state'] = 'None'
mdict = {'save_coefs': save_coefs, 'save_r': save_r, 'save_r2': save_r2, 'save_mse': save_mse,
         'save_iter': save_iter, 'totalTime': time.time() - time_start, 'nan_counter': nan_counter,
         'max_iter_counter': max_iter_counter, 'params': params}
if params['flag_tune'] > 0:
    mdict.update({'cv_results': cv_results, 'hyperparam_grid': hyperparam_grid})
if params['model_type'] in ['recon', 'decoding']:
    if params['flag_log'] == 1:
        save_y_pred = [np.exp(x.astype('float')) for x in save_y_pred]
    mdict.update({'save_y_pred': save_y_pred})
    if len(params['fixed_test']) == 0:
        mdict.update({'save_split_indices': save_split_indices, 'save_split_ratios': save_split_ratios,
                      'groups': groups})
sio.savemat(working_dir + fname_out, mdict=mdict)
if params['random_state'] == 'None':
    params['random_state'] = None
