%% initialization
if ~exist('DataB','var')
    load(data_file);
end

%% ---------------------------- THEIR CODE ----------------------------- %%
Features=NeuronType_Striatum(DataB);
Neurons=unique(DataB.Info.NeuronNumb,'rows');
AllNeurons = Neurons;
[Neurons,AllRamps,StereoCrit,MeanFR]=Selectramp(DataB,Neurons);
[~,AllStereo]=Selectstereo(DataB,Neurons);

%% epoch indices related to their output stucture
s1_idcs = 1 : 4;
s2_idcs = 5 : 8;	% what they say they do in the methods
% s2_idcs = 4 : 10;	% what's in their code
go_idcs = 9 : 10;

%% parse their ramp clusters

% preallocation
ramp_idcs = struct();
ramp_idcs_ud = struct();
nonramp_idcs = struct();
nonramp_idcs_ud = struct();

% S1-aligned "ramps"
% ramp_idcs.s1_onset.up = find(sum(AllRamps(:,1),2)>0);
% ramp_idcs.s1_onset.down = find(sum(AllRamps(:,2),2)>0);
% ramp_idcs.s1_onset.both = unique([...
%     ramp_idcs.s1_onset.up;...
%     ramp_idcs.s1_onset.down]);
% ramp_idcs.s1_offset.up = find(sum(AllRamps(:,3),2)>0);
% ramp_idcs.s1_offset.down = find(sum(AllRamps(:,4),2)>0);
% ramp_idcs.s1_offset.both = unique([...
%     ramp_idcs.s1_offset.up;...
%     ramp_idcs.s1_offset.down]);
% ramp_idcs.s1.up = unique([...
%     ramp_idcs.s1_onset.up;...
%     ramp_idcs.s1_offset.up]);
% ramp_idcs.s1.down = unique([...
%     ramp_idcs.s1_onset.down;...
%     ramp_idcs.s1_offset.down]);
% ramp_idcs.s1.both = unique([...
%     ramp_idcs.s1.up;...
%     ramp_idcs.s1.down]);
ramp_idcs_ud.S1_onset = find(sum(AllRamps(:,1:2),2)>0);
ramp_idcs_ud.S1_offset = find(sum(AllRamps(:,3:4),2)>0);
nonramp_idcs_ud.S1_onset = find(sum(AllRamps(:,1:2),2)==0);
nonramp_idcs_ud.S1_offset = find(sum(AllRamps(:,3:4),2)==0);

ramp_idcs.s1 = find(sum(AllRamps(:,s1_idcs),2)>0);
nonramp_idcs.s1 = find(sum(AllRamps(:,s1_idcs),2)==0);

% S2-aligned "ramps"
% ramp_idcs.s2_onset.up = find(sum(AllRamps(:,5),2)>0);
% ramp_idcs.s2_onset.down = find(sum(AllRamps(:,6),2)>0);
% ramp_idcs.s2_onset.both = unique([...
%     ramp_idcs.s2_onset.up;...
%     ramp_idcs.s2_onset.down]);
% ramp_idcs.s2_offset.up = find(sum(AllRamps(:,7),2)>0);
% ramp_idcs.s2_offset.down = find(sum(AllRamps(:,8),2)>0);
% ramp_idcs.s2_offset.both = unique([...
%     ramp_idcs.s2_offset.up;...
%     ramp_idcs.s2_offset.down]);
% ramp_idcs.s2.up = unique([...
%     ramp_idcs.s2_onset.up;...
%     ramp_idcs.s2_offset.up]);
% ramp_idcs.s2.down = unique([...
%     ramp_idcs.s2_onset.down;...
%     ramp_idcs.s2_offset.down]);
% ramp_idcs.s2.both = unique([...
%     ramp_idcs.s2.up;...
%     ramp_idcs.s2.down]);
ramp_idcs_ud.S2_onset = find(sum(AllRamps(:,5:6),2)>0);
ramp_idcs_ud.S2_offset = find(sum(AllRamps(:,7:8),2)>0);
nonramp_idcs_ud.S2_onset = find(sum(AllRamps(:,5:6),2)==0);
nonramp_idcs_ud.S2_offset = find(sum(AllRamps(:,7:8),2)==0);

ramp_idcs.s2 = find(sum(AllRamps(:,s2_idcs),2)>0);
nonramp_idcs.s2 = find(sum(AllRamps(:,s2_idcs),2)==0);

% go-aligned "ramps"
% ramp_idcs.go_cue.up = find(sum(AllRamps(:,9),2)>0);
% ramp_idcs.go_cue.down = find(sum(AllRamps(:,10),2)>0);
% ramp_idcs.go_cue.both = unique([...
%     ramp_idcs.go_cue.up;...
%     ramp_idcs.go_cue.down]);
ramp_idcs_ud.Go = find(sum(AllRamps(:,9:10),2)>0);
nonramp_idcs_ud.Go = find(sum(AllRamps(:,9:10),2)==0);

ramp_idcs.go = find(sum(AllRamps(:,go_idcs),2)>0);
nonramp_idcs.go = find(sum(AllRamps(:,go_idcs),2)==0);

% update "ramping" neurons
if exist('flagged_neurons','var')
    epochfields = fieldnames(ramp_idcs);
    n_epochfields = numel(epochfields);
    for ii = 1 : n_epochfields
        ramp_flags = ...
            ismember(ramp_idcs.(epochfields{ii}),flagged_neurons);
        ramp_idcs.(epochfields{ii}) = ...
            ramp_idcs.(epochfields{ii})(ramp_flags);
        nonramp_flags = ...
            ismember(nonramp_idcs.(epochfields{ii}),flagged_neurons);
        nonramp_idcs.(epochfields{ii}) = ...
            nonramp_idcs.(epochfields{ii})(nonramp_flags);
    end
    epochfields = fieldnames(ramp_idcs_ud);
    n_epochfields = numel(epochfields);
    for ii = 1 : n_epochfields
        ramp_flags = ...
            ismember(ramp_idcs_ud.(epochfields{ii}),flagged_neurons);
        ramp_idcs_ud.(epochfields{ii}) = ...
            ramp_idcs_ud.(epochfields{ii})(ramp_flags);
        nonramp_flags = ...
            ismember(nonramp_idcs_ud.(epochfields{ii}),flagged_neurons);
        nonramp_idcs_ud.(epochfields{ii}) = ...
            nonramp_idcs_ud.(epochfields{ii})(nonramp_flags);
    end
end

%% parse their non-stereotypical clusters

% preallocation
stereo_idcs = struct();
stereo_idcs_ud = struct();
nonstereo_idcs = struct();
nonstereo_idcs_ud = struct();

% S1-aligned non-stereotypical neurons
% nonstereo_idcs.s1_onset.up = find(sum(AllStereo(:,1),2)>0);
% nonstereo_idcs.s1_onset.down = find(sum(AllStereo(:,2),2)>0);
% nonstereo_idcs.s1_onset.both = unique([...
%     nonstereo_idcs.s1_onset.up;...
%     nonstereo_idcs.s1_onset.down]);
% nonstereo_idcs.s1_offset.up = find(sum(AllStereo(:,3),2)>0);
% nonstereo_idcs.s1_offset.down = find(sum(AllStereo(:,4),2)>0);
% nonstereo_idcs.s1_offset.both = unique([...
%     nonstereo_idcs.s1_offset.up;...
%     nonstereo_idcs.s1_offset.down]);
% nonstereo_idcs.s1.up = unique([...
%     nonstereo_idcs.s1_onset.up;...
%     nonstereo_idcs.s1_offset.up]);
% nonstereo_idcs.s1.down = unique([...
%     nonstereo_idcs.s1_onset.down;...
%     nonstereo_idcs.s1_offset.down]);
% nonstereo_idcs.s1.both = unique([...
%     nonstereo_idcs.s1.up;...
%     nonstereo_idcs.s1.down]);
stereo_idcs_ud.S1_onset = find(sum(AllStereo(:,1:2),2)>0);
stereo_idcs_ud.S1_offset = find(sum(AllStereo(:,3:4),2)>0);
nonstereo_idcs_ud.S1_onset = find(sum(AllStereo(:,1:2),2)==0);
nonstereo_idcs_ud.S1_offset = find(sum(AllStereo(:,3:4),2)==0);

stereo_idcs.s1 = find(sum(AllStereo(:,s1_idcs),2)>0);
nonstereo_idcs.s1 = find(sum(AllStereo(:,s1_idcs),2)==0);

% S2-aligned non-stereotypical neurons
% nonstereo_idcs.s2_onset.up = find(sum(AllStereo(:,5),2)>0);
% nonstereo_idcs.s2_onset.down = find(sum(AllStereo(:,6),2)>0);
% nonstereo_idcs.s2_onset.both = unique([...
%     nonstereo_idcs.s2_onset.up;...
%     nonstereo_idcs.s2_onset.down]);
% nonstereo_idcs.s2_offset.up = find(sum(AllStereo(:,7),2)>0);
% nonstereo_idcs.s2_offset.down = find(sum(AllStereo(:,8),2)>0);
% nonstereo_idcs.s2_offset.both = unique([...
%     nonstereo_idcs.s2_offset.up;...
%     nonstereo_idcs.s2_offset.down]);
% nonstereo_idcs.s2.up = unique([...
%     nonstereo_idcs.s2_onset.up;...
%     nonstereo_idcs.s2_offset.up]);
% nonstereo_idcs.s2.down = unique([...
%     nonstereo_idcs.s2_onset.down;...
%     nonstereo_idcs.s2_offset.down]);
% nonstereo_idcs.s2.both = unique([...
%     nonstereo_idcs.s2.up;...
%     nonstereo_idcs.s2.down]);
stereo_idcs_ud.S2_onset = find(sum(AllStereo(:,5:6),2)>0);
stereo_idcs_ud.S2_offset = find(sum(AllStereo(:,7:8),2)>0);
nonstereo_idcs_ud.S2_onset = find(sum(AllStereo(:,5:6),2)==0);
nonstereo_idcs_ud.S2_offset = find(sum(AllStereo(:,7:8),2)==0);

stereo_idcs.s2 = find(sum(AllStereo(:,s2_idcs),2)>0);
nonstereo_idcs.s2 = find(sum(AllStereo(:,s2_idcs),2)==0);

% go-aligned non-stereotypical neurons
% nonstereo_idcs.go_cue.up = find(sum(AllStereo(:,9),2)>0);
% nonstereo_idcs.go_cue.down = find(sum(AllStereo(:,10),2)>0);
% nonstereo_idcs.go_cue.both = unique([...
%     nonstereo_idcs.go_cue.up;...
%     nonstereo_idcs.go_cue.down]);
stereo_idcs_ud.Go = find(sum(AllStereo(:,9:10),2)>0);
nonstereo_idcs_ud.Go = find(sum(AllStereo(:,9:10),2)==0);

stereo_idcs.go = find(sum(AllStereo(:,go_idcs),2)>0);
nonstereo_idcs.go = find(sum(AllStereo(:,go_idcs),2)==0);

% update stereotypical neurons
if exist('flagged_neurons','var')
    epochfields = fieldnames(stereo_idcs);
    n_epochfields = numel(epochfields);
    for ii = 1 : n_epochfields
        flags = ismember(stereo_idcs.(epochfields{ii}),flagged_neurons);
        stereo_idcs.(epochfields{ii}) = stereo_idcs.(epochfields{ii})(flags);
    end
    epochfields = fieldnames(nonstereo_idcs);
    n_epochfields = numel(epochfields);
    for ii = 1 : n_epochfields
        flags = ismember(nonstereo_idcs.(epochfields{ii}),flagged_neurons);
        nonstereo_idcs.(epochfields{ii}) = nonstereo_idcs.(epochfields{ii})(flags);
    end
end

%% sanity check
A = [numel(ramp_idcs.s1),numel(nonramp_idcs.s1);...
numel(ramp_idcs.s2),numel(nonramp_idcs.s2);...
numel(ramp_idcs.go),numel(nonramp_idcs.go);...
numel(stereo_idcs.s1),numel(nonstereo_idcs.s1);...
numel(stereo_idcs.s2),numel(nonstereo_idcs.s2)];
[A,sum(A,2),A(:,1)./A(:,2)]