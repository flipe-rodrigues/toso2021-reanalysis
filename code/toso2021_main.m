%% initialization
warning('off');
close all;
clear;
clc;

%% task selection
task_str = 'duration';
% task_str = 'intensity';

%% directory settings
code_path = fileparts(matlab.desktop.editor.getActiveFilename);
addpath(genpath(code_path));
cd(code_path);
root_path = fileparts(code_path);
data_path = fullfile(root_path,'data');
data_filename = sprintf('%c%s_rats_subjectID.mat',...
    upper(task_str(1)),task_str(2:end));
data_file = fullfile(data_path,data_filename);
load(data_file);
data = DataB.Info;
clear DataB;

%% load reaction time data (not recorded simultaneously with DLS data)
rt_file = fullfile(data_path,'RT_behavior.mat');
rts = load(rt_file);
data.Rts = rts.(sprintf('%s_rats',capitalize(task_str)));

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

%% compute spike density functions
toso2021_computeSDFs;

%% neuron selection
toso2021_neuronSelection;

%% contrast settings

% contrast selection
% contrast_str = 't1';
% contrast_str = 't2';
% contrast_str = 'i1';
contrast_str = 'i2';
% contrast_str = 'choice';
% contrast_str = 'correct';

% parse selected contrast
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str,'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);
contrast_units = eval([contrast_str,'_units']);
contrast_lbl = upper(contrast_str);

%% figures 1-2
toso2021_eventDiagram;
toso2021_generalizationMatrix_Si;
toso2021_generalizationMatrix_Di;
toso2021_performance_Si;
toso2021_performance_Di;
toso2021_choiceGLM;
toso2021_psychometricCurves;

%% figures 3-5
toso2021_averageActivity_s2Aligned;
toso2021_rasters_s2Aligned;
toso2021_PCA_s2Aligned;
toso2021_neurometricCurves_s2Aligned;
toso2021_naiveBayesDecoder_avg_s2Aligned;

%% figure 6
toso2021_trialTypeDistributions;
toso2021_spikeCountGLM_data;
toso2021_positiveControlResponses;
toso2021_spikeCountGLM_simulations;

%% figure 7
toso2021_rampClustering;
toso2021_rampProportions;
toso2021_rampTemporalTuning;
toso2021_rampFiringRateRange;
toso2021_rampStereotypyCoefficients;
toso2021_rampBumpModelSpecification;
toso2021_rampDecodingPerformance_simulations;
toso2021_rampDecodingPerformance_s1s2;
toso2021_rampClusterability_data;
toso2021_rampClusterability_hybrid;

%% figure S1
toso2021_samplingScheme;

%% figure S2
toso2021_idealizedGeneralizationMatrix_Si;
toso2021_idealizedGeneralizationMatrix_Di;
toso2021_idealizedPerformance_Si;
toso2021_idealizedPerformance_Di;
toso2021_idealizedChoiceGLM;

%% figure S3
toso2021_psychometricModel;
toso2021_modelComparison;

%% figure S4
toso2021_averageActivity_s1Aligned;
toso2021_rasters_s1Aligned;
toso2021_PCA_s1Aligned;
toso2021_neurometricCurves_s1Aligned;
toso2021_naiveBayesDecoder_avg_s1Aligned;

%% figure S5
toso2021_correctnessIssues;