%% initialization
warning('off');
close all;
clear;
clc;

%% handle dependencies
% [file_list,~] = matlab.codetools.requiredFilesAndProducts(mfilename);
% n_files = numel(file_list);
% for ii = 1 : n_files
%
%     dest_path = fullfile(root_path,'dependencies');
%     copyfile(file_list{ii},dest_path);
% end

%% task selection
task_str = 'duration';
% task_str = 'intensity';

%% directory settings
root_path = fileparts(which(mfilename));
addpath(genpath(root_path));
cd(root_path);
data_path = fullfile(root_path,'data');
file_name = sprintf('%c%s_rats.mat',upper(task_str(1)),task_str(2:end));
data_file = fullfile(data_path,file_name);
% data_file = fullfile(data_path,'duration','AT4.mat');
load(data_file);
data = DataB.Info;
% data.Rts = DataB.RT;
% data.Action = 1 - data.Action;
clear DataB;

%% save settings
panel_path = fullfile(root_path,'panels',task_str);
raster_path = fullfile(root_path,'rasters',task_str);
if ~exist(panel_path,'dir')
    mkdir(panel_path);
end
if ~exist(raster_path,'dir')
    mkdir(raster_path);
end
want2save = true;

%% preface
toso2021_preface;

%% contrast settings
contrast_str = 't2';
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str(1:end-1),'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);
contrast_units = eval([contrast_str(1),'_units']);
contrast_lbl = [upper(contrast_str(1)),'_',contrast_str(2)];

contrasts = prev_t2;
contrast_set = unique(contrasts(valid_flags & ~isnan(contrasts)));
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = spring(n_contrasts);
contrast_lbl = ['prev. ',upper(contrast_str(1)),'_',contrast_str(2)];

%% script execution order
toso2021_contractionBias;
toso2021_choiceGLM;
toso2021_GLM_choice_categorical;
toso2021_psychometricCurves; return;
toso2021_trialTypeDistributions;
toso2021_neuronSelection;
toso2021_overallModulation;
toso2021_tiling;
toso2021_t2AlignedPCA;
toso2021_PCA;
toso2021_hierarchicalClustering;
toso2021_neurometricCurves;
toso2021_naiveBayesDecoder;
% toso2021_rasters;
return;

%% multiplicity
mult_flags = valid_flags;
stimuli = t2 * w2 + t1 * w1;
thresh = psycurves(contrast_mode_idx).fit.Fit(1) * range(stimuli) + min(stimuli)
t2_minus_t1 = stimuli(mult_flags) > thresh;
t2_plus_t1 = t2(mult_flags) + t1(mult_flags) > ...
    nanmedian(unique((t2(mult_flags) + t1(mult_flags))));
t2_alone = t2(mult_flags) > nanmedian(unique((t2(mult_flags))));

ground_truth = t2(mult_flags) - t1(mult_flags) > 0;

figure; hold on;
win = 50;
kernel = expkernel('mus',win,'binwidth',1);
plot(conv(choices(mult_flags) == ground_truth,kernel.pdf,'valid'),...
    'linewidth',1.5);
plot(conv(t2_minus_t1 == ground_truth,kernel.pdf,'valid'));
plot(conv(t2_plus_t1 == ground_truth,kernel.pdf,'valid'));
plot(conv(t2_alone == ground_truth,kernel.pdf,'valid'));

ylabel('Proportion')
xlabel('Trial #')
xlim([win,10e3])
ylim([0,1])

legend({...
    'choices == T_2 - T_1',...
    'w_2 \times T_2 + w_1 \times T_1 == T_2 - T_1',...
    'T_2 + T_1 == T_2 - T_1',...
    'T_2 > med(T_2) == T_2 - T_1'})
