%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% data curation
fprintf('\nHETEROGENEITY FIXES:\n');

% intensity
fprintf('INTENSITY:\n');
fprintf('- before: (intensity | I1 count | I2 count)\n');
summary(categorical([data.Intensity1,data.Intensity2]))
data.Intensity1 = round(data.Intensity1);
data.Intensity2 = round(data.Intensity2);
data.Intensity1(...
    data.Intensity1 == 58 | ...
    data.Intensity1 == 108 | ...
    data.Intensity1 == 200) = nan;
data.Intensity2(...
    data.Intensity2 == 58 | ...
    data.Intensity2 == 108 | ...
    data.Intensity2 == 200) = nan;
data.Intensity1(data.Intensity1 == 79) = 80;
data.Intensity2(data.Intensity2 == 79) = 80;
data.Intensity1(data.Intensity1 == 149) = 147;
data.Intensity2(data.Intensity2 == 149) = 147;
fprintf('- after: (intensity | I1 count | I2 count)\n');
summary(categorical([data.Intensity1,data.Intensity2]))

% duration
fprintf('DURATION:\n');
fprintf('- before: (duration | T1 count | T2 count)\n');
summary(categorical([data.Duration1,data.Duration2]))
data.Duration1(...
    data.Duration1 == 205 | ...
    data.Duration1 == 264 | ...
    data.Duration1 == 423 | ...
    data.Duration1 == 545) = nan;
data.Duration2(...
    data.Duration2 == 205 | ...
    data.Duration2 == 264 | ...
    data.Duration2 == 423 | ...
    data.Duration2 == 545) = nan;
fprintf('- after: (duration | T1 count | T2 count)\n');
summary(categorical([data.Duration1,data.Duration2]))

%% hard-coded task parameters
inferred_t1t2_bug = 250;
pre_init_padding = 1e3 + inferred_t1t2_bug * 0;
inter_t1t2_delay = 2e3 + inferred_t1t2_bug * 1;
post_t2_delay = 500;

%% neuron selection criteria
n_trial_cutoff = 3;
mean_fr_cutoff = 1;

%% parse meta data (bhv)
t1 = data.Duration1;
t2 = data.Duration2;
t1_set = unique(t1);
t2_set = unique(t2);
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
choices = data.Action;
choice_set = unique(choices);
n_choices = numel(choice_set);
prev_choices = [nan;choices(1:end)];
correct = choices == (t2 > t1);
pre_t1_delay = data.PreDelay;
trial_idcs = data.Trial;

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
    pre_t1_delay == 500 & ...
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

%% color scheme
t1_clrs = cool(n_t);
t2_clrs = spring(n_t);
i1_clrs = winter(n_i);
i2_clrs = copper(n_i);
choices_clrs = [.1,.5,1; .85,.1,.2];