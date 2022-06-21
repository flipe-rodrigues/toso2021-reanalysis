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
load(data_file);
data = DataB.Info;
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

%% script execution order
toso2021_samplingScheme;
toso2021_choiceGLM;
toso2021_psychometricCurves;
toso2021_trialTypeDistributions;
toso2021_neuronSelection;
toso2021_overallModulation;
toso2021_PCA;
% toso2021_rasters;
% toso2021_neurometricCurves;
% toso2021_naiveBayesDecoder;