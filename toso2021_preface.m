%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% data curation
fprintf('\nHETEROGENEITY FIXES:\n');

%% intensity heterogeneity fixes
fprintf('INTENSITY:\n');
fprintf('- before: (intensity | I1 count | I2 count)\n');
summary(categorical([data.Intensity1,data.Intensity2]))

% round
data.Intensity1 = round(data.Intensity1);
data.Intensity2 = round(data.Intensity2);

% group close neighbors
data.Intensity1(data.Intensity1 == 79) = 80;
data.Intensity2(data.Intensity2 == 79) = 80;
data.Intensity1(data.Intensity1 == 149) = 147;
data.Intensity2(data.Intensity2 == 149) = 147;

% remove spurious i1 trial types
i1_set_bugged = unique(data.Intensity1);
for ii = 1 : numel(i1_set_bugged)
    i1_flags = data.Intensity1 == i1_set_bugged(ii);
    if sum(i1_flags) < 100
        data.Intensity1(i1_flags) = nan;
    end
end

% remove spurious i2 trial types
i2_set_bugged = unique(data.Intensity2);
for ii = 1 : numel(i2_set_bugged)
    i2_flags = data.Intensity2 == i2_set_bugged(ii);
    if sum(i2_flags) < 100
        data.Intensity2(i2_flags) = nan;
    end
end

fprintf('- after: (intensity | I1 count | I2 count)\n');
summary(categorical([data.Intensity1,data.Intensity2]))

%% duration heterogeneity fixes
fprintf('DURATION:\n');
fprintf('- before: (duration | T1 count | T2 count)\n');
summary(categorical([data.Duration1,data.Duration2]))

% round
data.Duration1 = round(data.Duration1);
data.Duration2 = round(data.Duration2);

% remove spurious t1 trial types
t1_set_bugged = unique(data.Duration1);
for ii = 1 : numel(t1_set_bugged)
    t1_flags = data.Duration1 == t1_set_bugged(ii);
    if sum(t1_flags) < 100
        data.Duration1(t1_flags) = nan;
    end
end

% remove spurious t2 trial types
t2_set_bugged = unique(data.Duration2);
for ii = 1 : numel(t2_set_bugged)
    t2_flags = data.Duration2 == t2_set_bugged(ii);
    if sum(t2_flags) < 100
        data.Duration2(t2_flags) = nan;
    end
end

fprintf('- after: (duration | T1 count | T2 count)\n');
summary(categorical([data.Duration1,data.Duration2]))

%% hard-coded task parameters
inferred_misalignment = 200;
pre_init_padding = 1e3;
isi = 2e3;
post_t2_delay = 500;

%% neuron selection criteria
trial_count_cutoff = 2;
mean_fr_cutoff = 1;
stability_cutoff = .5;

%% parse meta data (bhv)
t1 = data.Duration1;
t2 = data.Duration2;
t1_set = unique(t1(~isnan(t1)));
t2_set = unique(t2(~isnan(t2)));
prev_t1_set = t1_set;
prev_t2_set = t2_set;
n_t1 = numel(t1_set);
n_t2 = numel(t2_set);
i1 = data.Intensity1;
i2 = data.Intensity2;
i1_set = unique(i1(~isnan(i1)));
i2_set = unique(i2(~isnan(i2)));
n_i1 = numel(i1_set);
n_i2 = numel(i2_set);
t_set = intersect(t1_set,t2_set);
i_set = intersect(i1_set,i2_set);
t1_mode_idx = find(t_set == mode(t1));
t2_mode_idx = find(t_set == mode(t2));
i1_mode_idx = find(i_set == mode(i1));
i2_mode_idx = find(i_set == mode(i2));
i1_max_idx = find(i_set == max(i1));
i2_max_idx = find(i_set == max(i2));
n_t = numel(t_set);
n_i = numel(i_set);
choice = data.Action;
choice_set = unique(choice);
n_choices = numel(choice_set);
pre_t1_delay = data.PreDelay + inferred_misalignment;
trial_idcs = data.Trial;
subject_ids = data.Subject;
subject_set = unique(subject_ids(~isnan(subject_ids)));

% previous trial
prev_t1 = [nan;t1(1:end-1)];
prev_t2 = [nan;t2(1:end-1)];
prev_i1 = [nan;i1(1:end-1)];
prev_i2 = [nan;i2(1:end-1)];
prev_choices = [nan;choice(1:end-1)];

%% units
t1_units = 'ms';
t2_units = 'ms';
prev_t1_units = 'ms';
prev_t2_units = 'ms';
i1_units = 'mm.s^{-1}';
i2_units = 'mm.s^{-1}';
prev_i1_units = 'mm.s^{-1}';
prev_i2_units = 'mm.s^{-1}';
choice_units = 'a.u.';

%% color scheme
t1_clrs = cool(n_t);
t2_clrs = colorlerp([.25,.5,1; [1,1,1]*.25; [1,1,0]],n_t);
prev_t1_clrs = autumn(n_t);
prev_t2_clrs = spring(n_t);
i1_clrs = winter(n_i);
i2_clrs = copper(n_i);
choice_clrs = [.1,.5,1; .85,.1,.2];

%% task variant adaptations

% delayed duration comparison 
if strcmpi(task_str,'duration')
    
    % stimuli
    s1 = t1;
    s2 = t2;
    s_set = t_set;
    
    % distractors
    d1 = i1;
    d2 = i2;
    d_set = i_set;
    
    % labels
    s1_lbl = 'T_1';
    s2_lbl = 'T_2';
    d1_lbl = 'I_1';
    d2_lbl = 'I_2';
    s_units = t1_units;
    d_units = i1_units;
    
    % colors
    s1_clrs = t1_clrs;
    s2_clrs = t2_clrs;
    
% delayed intensity comparison    
elseif strcmpi(task_str,'intensity')
    
    % stimuli
    s1 = i1;
    s2 = i2;
    s_set = i_set;
    
    % distractors
    d1 = t1;
    d2 = t2;
    d_set = t_set;
    
    % labels
    s1_lbl = 'I_1';
    s2_lbl = 'I_2';
    d1_lbl = 'T_1';
    d2_lbl = 'T_2';
    s_units = i1_units;
    d_units = t1_units;
    
    % colors
    s1_clrs = i1_clrs;
    s2_clrs = i2_clrs;
end

% modality agnostic set mode indices
s1_mode_idx = find(s_set == mode(s1));
s2_mode_idx = find(s_set == mode(s2));
d1_mode_idx = find(d_set == mode(d1));
d2_mode_idx = find(d_set == mode(d2));

% correctness
correct = choice == (s2 > s1);

%% kernel settings
psthbin = 1;
kernel = gammakernel('peakx',50,'binwidth',psthbin);
n_paddedtimebins = size(data.FR,2);
n_timebins = n_paddedtimebins - kernel.nbins + 1;
n_tbins = max(t_set) * psthbin;
padded_time = ...
    (1 : psthbin : n_paddedtimebins * psthbin) - psthbin;
validtime_flags = ...
    padded_time >= padded_time(1) - kernel.paddx(1) & ...
    padded_time <= padded_time(end) - kernel.paddx(end) + psthbin;
valid_time = padded_time(validtime_flags);

%% parse meta data (ephys)
neuron_idcs = unique(data.NeuronNumb);
n_neurons_total = numel(neuron_idcs);

%% trial pre-selection
valid_flags = ...
    pre_t1_delay == 500 + inferred_misalignment & ...
    ismember(i1,i_set) & ...
    ismember(t1,t_set) & ...
    ismember(i2,i_set) & ...
    ismember(t2,t_set);

%% figure options
figopt.color = 'w';
figopt.inverthardcopy = 'off';
figopt.numbertitle = 'off';
figopt.units = 'pixels';
figopt.menubar = 'figure';
figopt.toolbar = 'auto';

%% axes options

% default axes properties
axesopt.default.plotboxaspectratio = [1,1,1];
axesopt.default.ticklength = [1,1] * .025;
axesopt.default.linewidth = 2;
axesopt.default.xlimmode = 'auto';
axesopt.default.ylimmode = 'auto';
axesopt.default.xticklabelrotation = 0;
axesopt.default.yticklabelrotation = 0;
axesopt.default.fontname = 'helvetica';
axesopt.default.fontsize = 12;
axesopt.default.color = 'none';
axesopt.default.xcolor = 'k';
axesopt.default.ycolor = 'k';
axesopt.default.nextplot = 'add';
axesopt.default.tickdir = 'out';
axesopt.default.box = 'off';
axesopt.default.layer = 'top';

% psychometric-specific axes properties
axesopt.psycurve.ylim = [-.05,1.05];
axesopt.psycurve.ytick = linspace(0,1,3);
axesopt.psycurve.yticklabel = num2cell(axesopt.psycurve.ytick);
axesopt.psycurve.yticklabel(2) = {''};

% south-east inset properties (southeast)
axesopt.inset.se.position = [0.13+0.775*4/7,0.25,0.775/3,0.815/3];
axesopt.inset.se.fontsize = axesopt.default.fontsize * 5/6;
axesopt.inset.se.xaxislocation = 'bottom';
axesopt.inset.se.yaxislocation = 'right';

% color bar properties
axesopt.colorbar.ticklength = .025;
axesopt.colorbar.linewidth = 2;
axesopt.colorbar.fontname = 'helvetica';
axesopt.colorbar.fontsize = 12;
axesopt.colorbar.color = 'none';
axesopt.colorbar.xcolor = 'k';
axesopt.colorbar.ycolor = 'k';
axesopt.colorbar.tickdir = 'out';
axesopt.colorbar.box = 'off';