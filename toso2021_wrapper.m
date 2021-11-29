%% initialization
clear
clc;

%% task selection
task = 'duration';

%% directory settings
root_path = fileparts(which(mfilename));
addpath(genpath(root_path));
cd(root_path);
data_path = fullfile(root_path,'data');
file_name = sprintf('%c%s_rats.mat',upper(task(1)),task(2:end));
data_file = fullfile(data_path,file_name);
load(data_file);
data = DataB.Info;
clear DataB;

%% save settings
panel_path = fullfile(root_path,'panels');
raster_path = fullfile(root_path,'rasters');
if ~exist(panel_path,'dir')
    mkdir(panel_path);
end
if ~exist(raster_path,'dir')
    mkdir(raster_path);
end
want2save = true;

raster_path = 'C:\Users\flipe\Desktop\Comment - Toso et al. (2021)\rasters';

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