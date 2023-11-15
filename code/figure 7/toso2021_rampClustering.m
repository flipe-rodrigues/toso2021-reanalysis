%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% USE THEIR CODE TO CLASSIFY NEURONS AS "RAMPS" & "NON-RAMPS"
if ~exist('AllRamps','var')
    load(data_file);
    [Neurons,AllRamps,~,~,preselected_idcs] = Selectramp_THEIRCODE(DataB);
    clear DataB;
end

%% cluster settings
cluster_labels = {'ramp','nonramp'};
n_clusters = numel(cluster_labels);

%% epoch indices related to their output stucture

% preallocation
epoch_idcs = struct();

% epoch indices
epoch_idcs.s1_onset = 1 : 2;
epoch_idcs.s1_offset = 3 : 4;
epoch_idcs.s2_onset = 5 : 6;
epoch_idcs.s2_offset = 7 : 8;
epoch_idcs.go_cue = 9 : 10;
epoch_idcs.s1 = [epoch_idcs.s1_onset,epoch_idcs.s1_offset];
epoch_idcs.s2 = [epoch_idcs.s2_onset,epoch_idcs.s2_offset];

%% cluster epoch settings
cluster_epochs = fieldnames(epoch_idcs);
n_cluster_epochs = numel(cluster_epochs);

%% cluster roi settings
cluster_roi = [-1,1] * 500;
cluster_roi_n_bins = range(cluster_roi);
cluster_roi_time = linspace(cluster_roi(1),cluster_roi(2),cluster_roi_n_bins);

%% parse their ramp clusters

% preallocation
cluster_idcs = cell(n_cluster_epochs,n_clusters);

% iterate through epochs
for ee = 1 : n_cluster_epochs
    epoch = cluster_epochs{ee};
    
    % iterate through clusters
    for kk = 1 : n_clusters
        cluster = cluster_labels{kk};
        rule = sum(AllRamps(:,epoch_idcs.(epoch)),2) > 0 == strcmpi(cluster,'ramp');
        their_idcs = preselected_idcs(rule);
        intersection_flags = ismember(their_idcs,flagged_neurons);
        cluster_idcs{ee,kk} = their_idcs(intersection_flags);
    end
end

% table conversion
cluster_idcs = cell2table(cluster_idcs',...
    'variablenames',cluster_epochs,...
    'rownames',cluster_labels);

%% same but splitting between up- & down-ramping neurons

% preallocation
ramp_idcs = cell(n_cluster_epochs,n_clusters+1);

% iterate through epochs
for ee = 1 : n_cluster_epochs
    epoch = cluster_epochs{ee};
    
    % iterate through clusters
    for kk = 1 : n_clusters
        if ismember(epoch,{'s1','s2'})
            idcs = kk + [0,2];
        else
            idcs = kk;
        end
        
        % up- & down-ramping neurons
        rule = sum(AllRamps(:,epoch_idcs.(epoch)(idcs)),2) > 0;
        their_idcs = preselected_idcs(rule);
        intersection_flags = ismember(their_idcs,flagged_neurons);
        ramp_idcs{ee,kk} = their_idcs(intersection_flags);
    end
    
    % non-ramping neurons
    rule = sum(AllRamps(:,epoch_idcs.(epoch)),2) == 0;
    their_idcs = preselected_idcs(rule);
    intersection_flags = ismember(their_idcs,flagged_neurons);
    ramp_idcs{ee,end} = their_idcs(intersection_flags);
end

% table conversion
ramp_idcs = cell2table(ramp_idcs',...
    'variablenames',cluster_epochs,...
    'rownames',{'up','down','non'});