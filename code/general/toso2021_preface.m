%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% data curation
fprintf('\nHETEROGENEITY / ROUNDING ERROR FIXES:\n');

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

%% convert intensity from standard deviation to mean speed units
intensity_tfun = @(x) round(x .* sqrt(2 / pi));

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
post_s2_delay = 500;

%% neuron selection criteria
trial_count_cutoff = 1;
mean_fr_cutoff = 1;
stability_cutoff = .15;

%% parse meta data (bhv)
t1 = data.Duration1;
t2 = data.Duration2;
t1_set = unique(t1(~isnan(t1)));
t2_set = unique(t2(~isnan(t2)));
prev_t1_set = t1_set;
prev_t2_set = t2_set;
n_t1 = numel(t1_set);
n_t2 = numel(t2_set);
i1 = intensity_tfun(data.Intensity1);
i2 = intensity_tfun(data.Intensity2);
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
pre_s1_delay = data.PreDelay + inferred_misalignment;
trial_idcs = data.Trial;
n_total_trials = numel(trial_idcs);
if isfield(data,'Subject')
    subjects = data.Subject;
else
    subjects = ones(n_total_trials,1);
end
subject_set = unique(subjects(~isnan(subjects)));
n_subjects = numel(subject_set);

%% categorical T1
t1_cat = zeros(n_total_trials,1);
t1_cat(t1 > t_set(t1_mode_idx)) = +1;
t1_cat(t1 < t_set(t1_mode_idx)) = -1;
t1_cat = categorical(t1_cat,[-1,0,+1],...
    {'T_1 < 334 ms','T_1 = 334 ms','T_1 > 334 ms'});
t1_cat_set = unique(t1_cat(~isundefined(t1_cat)));
n_t1_cat = numel(t1_cat_set);

%% interaction terms
t1t2 = t1 .* t2;
t1t2_set = unique(t1t2(~isnan(t1t2)));
n_t1t2 = numel(t1t2_set);

%% previous trial
prev_t1 = [nan;t1(1:end-1)];
prev_t2 = [nan;t2(1:end-1)];
prev_i1 = [nan;i1(1:end-1)];
prev_i2 = [nan;i2(1:end-1)];
prev_choice = [nan;choice(1:end-1)];

%% units
t1_units = 'ms';
t2_units = 'ms';
t1t2_units = 'ms^{2}';
prev_t1_units = 'ms';
prev_t2_units = 'ms';
i1_units = 'mm.s^{-1}';
i2_units = 'mm.s^{-1}';
prev_i1_units = 'mm.s^{-1}';
prev_i2_units = 'mm.s^{-1}';
choice_units = 'a.u.';
correct_units = 'a.u.';
t1_cat_units = '';

%% color scheme
t1_clrs = cool(n_t) * .95;
t2_clrs = summer(n_t) * .95;
t1t2_clrs = parula(n_t1t2);
prev_t1_clrs = autumn(n_t);
prev_t2_clrs = spring(n_t);
i1_clrs = winter(n_i);
i2_clrs = copper(n_i);
choice_clrs = [.1,.5,1; .85,.1,.2];
prevchoice_clrs = choice_clrs / 2;
reward_clrs = [.25,.25,.25; .25,.9,.8];
prevreward_clrs = reward_clrs / 2;
subject_clr = [1,1,1] * .75;
stim_clrs = [1,1,1] .* [.75; 0];
t1_cat_clrs = cool(n_t1_cat);
ramp_clrs = [.85,.05,.25; .15,.15,.15];
rampud_clrs = [ramp_clrs(1,:);ramp_clrs(1,:);ramp_clrs(2,:)];
corr_clrs = [1,.25,0; .35,.35,.35];

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
    s1_lbl = 'T1';
    s2_lbl = 'T2';
    d1_lbl = 'I1';
    d2_lbl = 'I2';
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
    s1_lbl = 'I1';
    s2_lbl = 'I2';
    d1_lbl = 'T1';
    d2_lbl = 'T2';
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
correct_set = unique(correct);
correct_clrs = reward_clrs;
prev_correct = [nan;correct(1:end-1)];

%% choice & correctness intersection
choice_correct = (choice * 2 - 1) .* (correct + 1);
choice_correct = categorical(choice_correct,[-2,-1,1,2],{...
    'T2<T1_{correct}',...
    'T2<T1_{incorrect}',...
    'T2>T1_{incorrect}',...
    'T2>T1_{correct}'});
choice_correct_set = unique(choice_correct);
n_choice_correct = numel(choice_correct_set);
choice_correct_clrs = colorlerp(choice_clrs,n_choice_correct);
choice_correct_clrs = choice_correct_clrs([1,3,2,4],:);
choice_correct_units = 'a.u.';

%% parse meta data (ephys)
neuron_idcs = unique(data.NeuronNumb);
n_neurons_total = numel(neuron_idcs);

%% down-sample original spike counts
downsampling_factor = 1;
psthbin_src = 1;
psthbin = psthbin_src * downsampling_factor;
n_timebins_src = size(data.FR,2);
n_timebins = n_timebins_src / downsampling_factor;
time_src = 1 : psthbin_src : n_timebins_src * psthbin_src;
time = 1 : psthbin : n_timebins * psthbin;
if (psthbin_src ~= psthbin) && (size(data.FR,2) == n_timebins_src)
    data.FR = ...
        data.FR(:,1:downsampling_factor:end) + ...
        data.FR(:,downsampling_factor:downsampling_factor:end);
end

%% kernel settings
kernel_peak_time = 50;
gamma_kernel = gammakernel('peakx',kernel_peak_time,'binwidth',psthbin);
n_paddedtimebins = size(data.FR,2);
n_timebins = n_paddedtimebins - gamma_kernel.nbins + 1;
n_tbins = max(t_set) / psthbin;
padded_time = time;
validtime_flags = ...
    padded_time >= padded_time(1) - gamma_kernel.paddx(1) & ...
    padded_time <= padded_time(end) - gamma_kernel.paddx(end) + psthbin;
valid_time = padded_time(validtime_flags);

%% trial pre-selection
valid_flags = ...
    pre_s1_delay == 500 + inferred_misalignment & ...
    ismember(i1,i_set) & ...
    ismember(t1,t_set) & ...
    ismember(i2,i_set) & ...
    ismember(t2,t_set);

%% normalized stimulus dimensions

% duration
ntd = round((t2 - t1) ./ (t2 + t1),2);
ntd_set = unique(ntd(valid_flags));
n_ntd = numel(ntd_set);

% intensity
nid = round((i2 - i1) ./ (i2 + i1),1);
nid_set = unique(nid(valid_flags));
n_nid = numel(nid_set);
    
% delayed duration comparison
if strcmpi(task_str,'duration')
    
    % stimuli
    nsd = ntd;
    nsd_set = ntd_set;
    n_nsd = n_ntd;
    
    % distractors
    ndd = nid;
    ndd_set = nid_set;
    n_ndd = n_nid;
    
    % labels
    nsd_lbl = 'NTD';
    ndd_lbl = 'NID';
    
% delayed intensity comparison
elseif strcmpi(task_str,'intensity')

    % stimuli
    nsd = nid;
    nsd_set = nid_set;
    n_nsd = n_nid;
    
    % distractors
    ndd = ntd;
    ndd_set = ntd_set;
    n_ndd = n_ntd;
    
    % labels
    nsd_lbl = 'NID';
    ndd_lbl = 'NTD';
end

%% flag unique trials
pseudosession_transition_flags = [diff(data.Trial) ~= 1; true];
n_pseudosession_transitions = sum(pseudosession_transition_flags);
n_pseudosession_trialcounts = data.Trial(pseudosession_transition_flags);

% preallocation
unique_flags = false(size(pseudosession_transition_flags));
prev_session_rows = [];
prev_trial_count = 0;

% iterate through neuron transitions
for ii = 1 : n_pseudosession_transitions
    if ii > 1
        prev_trial_count = prev_trial_count + ...
            n_pseudosession_trialcounts(ii-1);
    end
    pseudo_session_rows = ...
        (1 : n_pseudosession_trialcounts(ii)) + prev_trial_count;
    if ii == 0
        unique_flags(pseudo_session_rows) = true;
    else
        if numel(pseudo_session_rows) ~= numel(prev_session_rows)
            unique_flags(pseudo_session_rows) = true;
        elseif any(t1(pseudo_session_rows) ~= t1(prev_session_rows)) && ...
                any(i1(pseudo_session_rows) ~= i1(prev_session_rows)) && ...
                any(t2(pseudo_session_rows) ~= t2(prev_session_rows)) && ...
                any(i2(pseudo_session_rows) ~= i2(prev_session_rows))
            unique_flags(pseudo_session_rows) = true;
        end
    end
    prev_session_rows = pseudo_session_rows;
end

%% fade settings
fadeifnoisy = true;
alphabounds_sem = [.05,.25];
alphabounds_mu = [.15,1];

%% figure options
figopt.color = 'w';
figopt.inverthardcopy = 'off';
figopt.numbertitle = 'off';
figopt.units = 'pixels';
figopt.menubar = 'figure';
figopt.toolbar = 'auto';

%% axes options
axesopt = struct();

% default axes properties
axesopt.default.plotboxaspectratio = [1,1,1];
axesopt.default.ticklength = [1,1] * .025;
axesopt.default.linewidth = 2;
axesopt.default.xlimmode = 'auto';
axesopt.default.ylimmode = 'auto';
axesopt.default.xlimspec = 'tight';
axesopt.default.ylimspec = 'tight';
axesopt.default.xticklabelrotation = 0;
axesopt.default.yticklabelrotation = 0;
axesopt.default.fontname = 'helvetica';
axesopt.default.fontsize = 12;
axesopt.default.color = 'none';
axesopt.default.xcolor = 'k';
axesopt.default.ycolor = 'k';
axesopt.default.zcolor = 'k';
axesopt.default.nextplot = 'add';
axesopt.default.xdir = 'normal';
axesopt.default.ydir = 'normal';
axesopt.default.tickdir = 'out';
axesopt.default.box = 'off';
axesopt.default.layer = 'top';

% psychometric-specific axes properties
axesopt.psycurve.ylim = [-.05,1.05];
axesopt.psycurve.ytick = linspace(0,1,5);
axesopt.psycurve.yticklabel = num2cell(axesopt.psycurve.ytick);
axesopt.psycurve.yticklabel(~ismember(axesopt.psycurve.ytick,[0,1])) = {''};

% south-east inset properties (southeast)
axesopt.inset.se.position = [0.13+0.775*4/7,0.25,0.775/3,0.815/3];
axesopt.inset.se.fontsize = axesopt.default.fontsize * 5/6;
axesopt.inset.se.xaxislocation = 'bottom';
axesopt.inset.se.yaxislocation = 'right';

% south-east inset properties (southeast)
axesopt.inset.ne.position = [0.13+0.775*4/7,0.65,0.775/3,0.815/3];
axesopt.inset.ne.fontsize = axesopt.default.fontsize * 5/6;
axesopt.inset.ne.xaxislocation = 'bottom';
axesopt.inset.ne.yaxislocation = 'right';

% north-west inset properties
axesopt.inset.nw.position = [.275,.65,0.775/3,0.815/3];
axesopt.inset.nw.fontsize = axesopt.default.fontsize * 5/6;
axesopt.inset.nw.xaxislocation = 'top';
axesopt.inset.nw.yaxislocation = 'left';

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