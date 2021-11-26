%% initialization
clear
clc;

%% task selection
task = 'duration';

%% data directory settings
file_name = sprintf('%c%s_rats.mat',upper(task(1)),task(2:end));
data_file = fullfile(pwd,file_name);
load(data_file);
data = DataB.Info;
clear DataB;

%% save settings
save_path = 'C:\Users\flipe\Desktop\Comment - Toso et al. (2021)';
mkdir(fullfile(save_path,'panels'));
mkdir(fullfile(save_path,'rasters'));
want2save = true;

%% script execution order
toso2021_preface;
toso2021_behavior;
toso2021_neuronSelection;
toso2021_trialTypeDistributions;
toso2021_overallModulation;
toso2021_PCA;
% toso2021_rasters;
toso2021_neurometricCurves;
% toso2021_naiveBayesDecoder;