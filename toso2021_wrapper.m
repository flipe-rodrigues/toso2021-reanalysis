%% initialization
clear
clc;

%% task selection
task = 'duration';

%% data directory settings
root_path = pwd;
data_path = fullfile(root_path,'data');
file_name = sprintf('%c%s_rats.mat',upper(task(1)),task(2:end));
data_file = fullfile(data_path,file_name);
load(data_file);
data = DataB.Info;
clear DataB;

%% save settings
panel_path = fullfile(root_path,'panels');
raster_path = fullfile(root_path,'rasters');
mkdir(panel_path);
mkdir(raster_path);
want2save = true;

%% script execution order
toso2021_preface;
toso2021_behavior;
toso2021_trialTypeDistributions;
toso2021_neuronSelection;
toso2021_overallModulation;
toso2021_PCA;
% toso2021_rasters;
toso2021_neurometricCurves;
% toso2021_naiveBayesDecoder;