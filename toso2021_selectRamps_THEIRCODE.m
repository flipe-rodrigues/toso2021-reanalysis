%% initialization
if ~exist('DataB','var')
    load(data_file);
end

%% ---------------------------- THEIR CODE ----------------------------- %%
Features = NeuronType_Striatum(DataB);
Neurons=unique(DataB.Info.NeuronNumb,'rows');
AllNeurons = Neurons;
[Neurons,AllRamps,~,~,preselected_idcs] = Selectramp(DataB,Neurons);
[Neurons,AllCorrs] = Selectcorr(DataB,Neurons);

%% epoch indices related to their output stucture
s1_idcs = 1 : 4;
s2_idcs = 5 : 8;	% what they say they do in the methods
% s2_idcs = 4 : 10;	% what's in their code
go_idcs = 9 : 10;

%% parse their ramp clusters

% preallocation
ramp_idcs = struct();
nonramp_idcs = struct();

% S1-aligned "ramps"
ramp_idcs.s1 = preselected_idcs(sum(AllRamps(:,s1_idcs),2)>0);
nonramp_idcs.s1 = preselected_idcs(sum(AllRamps(:,s1_idcs),2)==0);

% S2-aligned "ramps"
ramp_idcs.s2 = preselected_idcs(sum(AllRamps(:,s2_idcs),2)>0);
nonramp_idcs.s2 = preselected_idcs(sum(AllRamps(:,s2_idcs),2)==0);

% go-aligned "ramps"
ramp_idcs.go = preselected_idcs(sum(AllRamps(:,go_idcs),2)>0);
nonramp_idcs.go = preselected_idcs(sum(AllRamps(:,go_idcs),2)==0);

%% parse their correlated clusters

% preallocation
corr_idcs = struct();
noncorr_idcs = struct();

% S1-aligned "corrs"
corr_idcs.s1 = preselected_idcs(sum(AllCorrs(:,s1_idcs),2)>0);
noncorr_idcs.s1 = preselected_idcs(sum(AllCorrs(:,s1_idcs),2)==0);

% S2-aligned "corrs"
corr_idcs.s2 = preselected_idcs(sum(AllCorrs(:,s2_idcs),2)>0);
noncorr_idcs.s2 = preselected_idcs(sum(AllCorrs(:,s2_idcs),2)==0);

% go-aligned "corrs"
corr_idcs.go = preselected_idcs(sum(AllCorrs(:,go_idcs),2)>0);
noncorr_idcs.go = preselected_idcs(sum(AllCorrs(:,go_idcs),2)==0);