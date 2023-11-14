%% initialization
warning('off');
close all;
clear;
clc;

%% task selection
task_str = 'duration';
% task_str = 'intensity';

%% directory settings
code_path = fileparts(which(mfilename));
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
X = data.FR;
g = circshift(gamma_kernel.pdf,gamma_kernel.nbins/2+1);
g = padarray(g,[0,n_paddedtimebins-gamma_kernel.nbins],0,'post');
G = cell2mat(arrayfun(@(i)circshift(g,i),(1:n_paddedtimebins)'-1,...
    'uniformoutput',false));
Z = X * G / psthbin * 1e3;
data.SDF = Z;
clear X G Z;

%% neuron selection
toso2021_neuronSelection;
% toso2021_parseSessionData;

%% cluster neurons into ramping & non-ramping neurons using their code
toso2021_selectRamps;

%% contrast selection
contrast_str = 'i2';
% contrast_str = 't2';
% contrast_str = 'i1';
% contrast_str = 't1';
% contrast_str = 'choice';

%% parse selected contrast
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str,'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);
contrast_units = eval([contrast_str,'_units']);
contrast_lbl = [upper(contrast_str(1)),'_',contrast_str(2)];

%% figures 1 & 2
toso2021_eventDiagram;
toso2021_generalizationMatrix_Si;
toso2021_generalizationMatrix_Di;
toso2021_performance_Si;
toso2021_performance_Di;
toso2021_choiceGLM;
toso2021_psychometricCurves;

%% figures 3-5
toso2021_averageActivity_s2Aligned;
toso2021_rasters_s2Aligned_i1_t1_i2_choice;
toso2021_PCA_s2Aligned;
toso2021_neurometricCurves_s2Aligned;
toso2021_naiveBayesDecoder_avg_s2Aligned;

%% figure 6
% toso2021_spikeCountsGLM_data;
% toso2021_spikeCountsGLM_simulations;

%% figure 7


%% figure S1
toso2021_samplingScheme;

%% figure S2
toso2021_averageActivity_s1Aligned;
toso2021_rasters_s1Aligned_i1_t1_i2_choice;
toso2021_PCA_s1Aligned;
toso2021_neurometricCurves_s1Aligned;
toso2021_naiveBayesDecoder_avg_s1Aligned;