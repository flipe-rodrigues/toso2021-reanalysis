%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% stimulus settings
stimuli = round(s2 - s1);
stim_set = unique(stimuli(valid_flags));
stim2group_flags = abs(diff(stim_set)) <= range(stim_set) * .01;
if any(stim2group_flags)
    for ii = find(stim2group_flags)'
        stimuli(stimuli == stim_set(ii)) = stim_set(ii + 1);
    end
end
stim_set = unique(stimuli(valid_flags));
n_stimuli = numel(stim_set);
stim_lbl = sprintf('%s - %s (%s)',s2_lbl,s1_lbl,s_units);

%% construct psychophysical triple

% normalize stimulus range
norm_stimuli = (stimuli - min(stimuli)) / range(stimuli);
normstim_set = unique(norm_stimuli(valid_flags));

% preallocation
psycurves = struct();

% iterate through subjects
for ss = 1 : n_subjects
    subject_flags = subjects == subject_set(ss);
    
    % iterate through stimuli
    for ii = 1 : n_stimuli
        stim_flags = stimuli == stim_set(ii);
        trial_flags = ...
            valid_flags & ...
            subject_flags & ...
            stim_flags;
        psycurves(ss).x(ii,1) = normstim_set(ii);
        psycurves(ss).y(ii,1) = sum(choice(trial_flags));
        psycurves(ss).n(ii,1) = sum(trial_flags);
        psycurves(ss).err(ii,1) = ...
            std(choice(trial_flags)) / sqrt(sum(trial_flags));
    end
end

%% fit psychometric function

% psychometric fit settings
psyopt.fit = struct();
psyopt.fit.expType = 'YesNo';
psyopt.fit.sigmoidName = 'logistic';
psyopt.fit.estimateType = 'MAP';
psyopt.fit.confP = [.95,.9,.68];
psyopt.fit.borders = [0,1;0,1;0,.25;0,.25;0,0];
psyopt.fit.fixedPars = [nan,nan,nan,nan,0];
psyopt.fit.stepN = [100,100,20,20,20];

% iterate through subjects
for ss = 1 : n_subjects
    
    % fit psychometric curve
    psycurves(ss).fit = ...
        psignifit([psycurves(ss).x,psycurves(ss).y,psycurves(ss).n],psyopt.fit);
end

%% plot phychometric function

% stimulus-specific axes properties
axesopt.stimulus.xlim = ...
    ([normstim_set(1),normstim_set(end)] +  [-1,1] * .05);
axesopt.stimulus.xtick = normstim_set;
axesopt.stimulus.xticklabel = num2cell(round(stim_set,2));
ticks2delete = ...
    ~ismember(axesopt.stimulus.xtick,...
    [min(axesopt.stimulus.xtick),max(axesopt.stimulus.xtick)]);
axesopt.stimulus.xticklabel(ticks2delete) = {''};

% psychometric plot settings
psyopt.plot = struct();
psyopt.plot.linewidth = 1.5;
psyopt.plot.marker = 'o';
psyopt.plot.markersize = 8.5;
psyopt.plot.plotdata = true;
psyopt.plot.gradeclrs = false;
psyopt.plot.patchci = false;
psyopt.plot.normalizemarkersize = false;

% figure initialization
fig = figure(figopt,...
    'name',sprintf('performance_t2t1'));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
xlabel(stim_lbl);
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*median(normstim_set),ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through subjects
for ss = n_subjects : -1 : 1
    clr = [0,0,0] + [1,1,1] * .75 * (ss ~= 1);
    
    % plot psychometric curve
    psyopt.plot.datafaceclr = clr;
    psyopt.plot.dataedgeclr = clr;
    psyopt.plot.overallvisibility = 'off';
    psyopt.plot.normalizemarkersize = false;
    psyopt.plot.plotfit = false;
    p(ss) = plotpsy(psycurves(ss),psycurves(ss).fit,psyopt.plot);
    plot(psycurves(ss).x,...
        psycurves(ss).y./psycurves(ss).n,...
        'color',clr,...
        'linestyle','--',...
        'linewidth',1.5);
    plot(psycurves(ss).x(1:n_stimuli/2),...
        psycurves(ss).y(1:n_stimuli/2)./psycurves(ss).n(1:n_stimuli/2),...
        'color',clr,...
        'linestyle','-',...
        'linewidth',1.5);
    plot(psycurves(ss).x(n_stimuli/2+1:end),...
        psycurves(ss).y(n_stimuli/2+1:end)./psycurves(ss).n(n_stimuli/2+1:end),...
        'color',clr,...
        'linestyle','-',...
        'linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot performance as a function of NTD

% figure initialization
fig = figure(figopt,...
    'name',sprintf('performance_NTD_distros'));

% axes initialization
axes(...
    axesopt.default,...
    'plotboxaspectratio',[1,2,1],...
    'xlim',[min(ntd_set),max(ntd_set)]+[-1,1]*1/2*range(ntd_set),...
    'xtick',ntd_set,...
    'ylim',axesopt.psycurve.ylim,...
    'ytick',axesopt.psycurve.ytick,...
    'yticklabel',axesopt.psycurve.yticklabel);
xlabel('NTD');
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through subjects
for ss = n_subjects : -1 : 1
    clr = [0,0,0] + [1,1,1] * .75 * (ss ~= 1);
    offset = (ss - (n_subjects + 1) / 2) * .025 * range(xlim);
    
    % plot performance
    plot(ntd_set(1+(psycurves(ss).x>.5))+offset,...
        psycurves(ss).y./psycurves(ss).n,...
        'marker','o',...
        'markersize',8.5,...
        'markeredgecolor',clr,...
        'markerfacecolor',clr,...
        'linestyle','none',...
        'linewidth',1.5);
    plot(ntd_set+offset,...
        psycurves(ss).y(floor(n_stimuli/2)+[0,1])./...
        psycurves(ss).n(floor(n_stimuli/2)+[0,1]),...
        'color',clr,...
        'linestyle','--',...
        'linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot performance as a function of NTD

% figure initialization
fig = figure(figopt,...
    'name',sprintf('performance_NTD_averages'));

% axes initialization
axes(...
    axesopt.default,...
    'plotboxaspectratio',[1,2,1],...
    'xlim',[min(ntd_set),max(ntd_set)]+[-1,1]*1/2*range(ntd_set),...
    'xtick',ntd_set,...
    'ylim',axesopt.psycurve.ylim,...
    'ytick',axesopt.psycurve.ytick,...
    'yticklabel',axesopt.psycurve.yticklabel);
xlabel('NTD');
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through subjects
for ss = n_subjects : -1 : 1
    clr = [0,0,0] + [1,1,1] * .75 * (ss ~= 1);
    offset = (ss - (n_subjects + 1) / 2) * .025 * range(xlim);
    
    % compute average performance per NTD level
    ntd_perf = [...
        psycurves(ss).y(1:floor(n_stimuli/2))./psycurves(ss).n(1:floor(n_stimuli/2)),...
        psycurves(ss).y(ceil(n_stimuli/2)+1:end)./psycurves(ss).n(ceil(n_stimuli/2)+1:end)];
    ntd_mus = mean(ntd_perf,1);
    ntd_std = std(ntd_perf,0,1);
    
    % plot performance
    errorbar(ntd_set+offset,ntd_mus,ntd_std,...
        'color',clr,...
        'marker','s',...
        'markersize',12.5,...
        'markeredgecolor','w',...
        'markerfacecolor',clr,...
        'linestyle','-',...
        'linewidth',1.5,...
        'capsize',3);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end