%% initialization
warning('off');
close all;
clear;
clc;

%% task selection
task_str = 'duration';
% task_str = 'intensity';

%% directory settings
root_path = fileparts(which(mfilename));
addpath(genpath(root_path));
cd(root_path);
data_path = fullfile(root_path,'data');
file_name = sprintf('%c%s_rats_ok.mat',upper(task_str(1)),task_str(2:end));
data_file = fullfile(data_path,file_name);
load(data_file);
data = DataB.Info;
% ---------------------------- THEIR RAMPS ------------------------------ %
% Features=NeuronType_Striatum(DataB);
% Neurons=unique(DataB.Info.NeuronNumb,'rows');
% AllNeurons=Neurons;
% [Neurons,AllRamps]=Selectramp(DataB,Neurons);
% ramp_idcs.s1.on.up = find(sum(AllRamps(:,1),2)>0);
% ramp_idcs.s1.on.down = find(sum(AllRamps(:,2),2)>0);
% ramp_idcs.s1.off.up = find(sum(AllRamps(:,3),2)>0);
% ramp_idcs.s1.off.down = find(sum(AllRamps(:,4),2)>0);
% ramp_idcs.s2.on.up = find(sum(AllRamps(:,5),2)>0);
% ramp_idcs.s2.on.down = find(sum(AllRamps(:,6),2)>0);
% ramp_idcs.s2.off.up = find(sum(AllRamps(:,7),2)>0);
% ramp_idcs.s2.off.down = find(sum(AllRamps(:,8),2)>0);
% ----------------------------------------------------------------------- %
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
downsampling_factor = 1;
kernel_peak_time = 25;
toso2021_preface;

%%

% timestamped task events
task_event_times = cumsum([...
    repmat(pre_init_padding,n_total_trials,1),...   % initiation
    pre_s1_delay,...                                % T1 onset
    t1,...                                          % T1 offset
    repmat(isi,n_total_trials,1),...                % T2 onset
    t2,...                                          % T2 offset
    repmat(post_s2_delay,n_total_trials,1)],2);     % go cue

%
unique_task_event_combs = unique(task_event_times(valid_flags,:),'rows');
[n_event_combs,n_event_types] = size(unique_task_event_combs);    

% kernel definition
K = normpdf(padded_time,padded_time',28);
K = K ./ nansum(K);

% preallocation
data.SDF = nan(n_total_trials,n_paddedtimebins);

% iterate through trials
for ii = 1 : n_event_combs
    progressreport(ii,n_event_combs,'event-informed SDF estimation');
    
%     figure('windowstyle','docked');
    
    event_bounds = [0,unique_task_event_combs(ii,:),n_paddedtimebins];
    trial_flags = ...
        valid_flags & ...
        all(task_event_times == unique_task_event_combs(ii,:),2);

    % iterate through events
    for jj = 1 : n_event_types + 1
        time_flags = ...
            padded_time >= event_bounds(jj) & ...
            padded_time < event_bounds(jj + 1);
        
        %
        X = data.FR(trial_flags,time_flags);
        k = K(time_flags,time_flags);
        k = k ./ nansum(k);
        Z = X * k;
        
%         figure;
%         hold on;
%         plot(padded_time(time_flags),k)
%         [~,idx] = min(abs(padded_time(time_flags) - mean(padded_time(time_flags))));
%         plot(padded_time(time_flags),k(1,:),'k','linewidth',1.5)
%         plot(padded_time(time_flags),k(idx,:),'k','linewidth',1.5)
%         a=1
        
%         subplot(2,1,1);
%         hold on;
%         imagesc([event_bounds(jj),event_bounds(jj+1)],[1,n_total_trials],X)
%         axis tight;
%         subplot(2,1,2);
%         hold on;
%         imagesc([event_bounds(jj),event_bounds(jj+1)],[1,n_total_trials],Z)
%         axis tight;
%         ax1=subplot(2,1,1);
%         ax2=subplot(2,1,2);
%         linkaxes([ax1,ax2]);
        
        %
        data.SDF(trial_flags,time_flags) = Z;
    end
    
    %
    X = data.FR(trial_flags,:);
    gamma_kernel = gammakernel('peakx',50,'binwidth',psthbin);
    g = circshift(gamma_kernel.pdf,gamma_kernel.nbins/2+1);
    g = padarray(g,[0,n_paddedtimebins-gamma_kernel.nbins],0,'post');
    G = cell2mat(arrayfun(@(i)circshift(g,i),(1:n_paddedtimebins)'-1,...
        'uniformoutput',false));
    Z = X * G / psthbin * 1e3;
    data.SDF(trial_flags,:) = Z;
end

%% contrast settings
contrast_str = 'i1';
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str,'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);
contrast_units = eval([contrast_str,'_units']);
contrast_lbl = [upper(contrast_str(1)),'_',contrast_str(2)];

%% script execution order

% behavior
toso2021_eventDiagram;
toso2021_generalizationMatrix_Si;
toso2021_generalizationMatrix_Di;
toso2021_performance_Si;
toso2021_performance_Di;
toso2021_GLM_choice;
toso2021_psychometricCurves;

% ephys
toso2021_neuronSelection;
return;
toso2021_averageActivity
toso2021_PCA_t2Aligned;
toso2021_neurometricCurves;
toso2021_naiveBayesDecoder;