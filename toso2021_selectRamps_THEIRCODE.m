%% initialization
if ~exist('DataB','var')
    load(data_file);
end

%% ---------------------------- THEIR CODE ----------------------------- %%
Features = NeuronType_Striatum(DataB);
Neurons=unique(DataB.Info.NeuronNumb,'rows');
AllNeurons = Neurons;
[Neurons,AllRamps,~,~,preselected_idcs] = Selectramp(DataB,Neurons);

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
cluster_roi_n_bins = range(cluster_roi) / psthbin;
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

%%

% preallocation
ramp_idcs = cell(n_cluster_epochs - 2,n_clusters);

% iterate through epochs
for ee = 1 : n_cluster_epochs - 2
    epoch = cluster_epochs{ee};

    % iterate through clusters
    for kk = 1 : 2
        rule = AllRamps(:,epoch_idcs.(epoch)(kk)) > 0;
        their_idcs = preselected_idcs(rule);
        intersection_flags = ismember(their_idcs,flagged_neurons);
        ramp_idcs{ee,kk} = their_idcs(intersection_flags);
    end
end

% table conversion
ramp_idcs = cell2table(ramp_idcs',...
    'variablenames',cluster_epochs(1:end-2),...
    'rownames',{'up','down'});