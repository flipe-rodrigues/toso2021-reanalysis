%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% stimulus settings
w_norm = sum(abs(beta_s1) + abs(beta_s2));
w1 = beta_s1 / w_norm;
w2 = beta_s2 / w_norm;
stimuli = round(s2 * w2 + s1 * w1);
stimuli = s2;
stim_set = unique(stimuli(valid_flags));
stim2group_flags = abs(diff(stim_set)) <= stim_set(2:end) * .05;
if any(stim2group_flags)
    for ii = find(stim2group_flags)'
        stimuli(stimuli == stim_set(ii)) = stim_set(ii + 1);
    end
end
stim_set = unique(stimuli(valid_flags));
stim_mode_idx = find(stim_set == mode(stimuli));
n_stimuli = numel(stim_set);
stim_lbl = sprintf('%.2f \\times %s + %.2f \\times %s (%s)',...
    w2,s2_lbl,w1,s1_lbl,s_units);

%% construct psychophysical triples

% normalize stimulus range
norm_stimuli = (stimuli - min(stimuli)) / range(stimuli);
normstim_set = unique(norm_stimuli(valid_flags));

% preallocation
psycurves = struct();

% iterate through contrasts
for kk = 1 : n_contrasts
    contrast_flags = contrasts == contrast_set(kk);
    
    % iterate through stimuli
    for ii = 1 : n_stimuli
        stimulus_flags = stimuli == stim_set(ii);
        trial_flags = ...
            valid_flags & ...
            unique_flags & ...
            contrast_flags & ...
            stimulus_flags;
        psycurves(kk).x(ii,1) = normstim_set(ii);
        psycurves(kk).y(ii,1) = sum(choice(trial_flags));
        psycurves(kk).n(ii,1) = sum(trial_flags);
        psycurves(kk).err(ii,1) = ...
            std(choice(trial_flags)) / sqrt(sum(trial_flags));
    end
end

% iterate through stimuli
for ii = 1 : n_stimuli
    stimulus_flags = stimuli == stim_set(ii);
    trial_flags = ...
        valid_flags & ...
        stimulus_flags;
    bigpsy.x(ii,1) = normstim_set(ii);
    bigpsy.y(ii,1) = sum(choice(trial_flags));
    bigpsy.n(ii,1) = sum(trial_flags);
    bigpsy.err(ii,1) = ...
        std(choice(trial_flags)) / sqrt(sum(trial_flags));
end

%% fit psychometric functions

% psychometric fit settings
psyopt.fit = struct();
psyopt.fit.expType = 'YesNo';
psyopt.fit.sigmoidName = 'rgumbel';
psyopt.fit.estimateType = 'MAP';
psyopt.fit.confP = [.95,.9,.68];
psyopt.fit.borders = [0,1; 0,1; 0,.25; 0,.25; 0,0];
psyopt.fit.fixedPars = [nan,nan,nan,nan,0];
psyopt.fit.stepN = [100,100,40,40,20];

% iterate through contrasts
for kk = 1 : n_contrasts

    % fit psychometric curve
    psycurves(kk).fit = ...
        psignifit([psycurves(kk).x,psycurves(kk).y,psycurves(kk).n],psyopt.fit);
end

%% plot phychometric functions

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
psyopt.plot.normalizemarkersize = true;

% figure initialization
fig = figure(figopt,...
    'name',sprintf('psychometric_curves_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
xlabel(sprintf('%s (%s)',s2_lbl,s_units));
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*median(normstim_set),ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through contrasts
for kk = 1 : n_contrasts

    % plot psychometric curve
    psyopt.plot.datafaceclr = contrast_clrs(kk,:);
    psyopt.plot.overallvisibility = 'off';
    psyopt.plot.plotfit = sum(psycurves(kk).n ~= 0) >= 2;
    p(kk) = plotpsy(psycurves(kk),psycurves(kk).fit,psyopt.plot);
end

% legend
if iscategorical(contrasts)
    leg_str = categories(contrast_set);
else
    leg_str = cellfun(@(x,y)sprintf('%s = %i %s',x,y,contrast_units),...
        repmat({contrast_lbl},n_contrasts,1),num2cell(contrast_set),...
        'uniformoutput',false);
end
legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
    'edgecolor','k',...
    'position',[0.085,0.66,.27,.2],...
    'box','on');

%% inset with delta P(T2 > T1)
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.inset.se,...
    'ylim',[-.25,.25]+[-1,1]*.05*.5,...
    'ytick',[-.25,0,.25]);
ylabel(sprintf('\\DeltaP(%s > %s)',s2_lbl,s1_lbl),...
    'rotation',-90,...
    'verticalalignment','bottom');

% categorical boundary
plot([1,1]*median(normstim_set),ylim,...
    'color','k',...
    'linestyle',':');

% zero level
plot(xlim,[0,0],...
    'color','k',...
    'linestyle',':');

% iterate through contrasts
for kk = 1 : n_contrasts
    if kk == contrast_mode_idx
        continue;
    end
    
    % plot delta P(T2 > T1)
    p_ctrl = psycurves(contrast_mode_idx).y ./ psycurves(contrast_mode_idx).n;
    p_cond = psycurves(kk).y ./ psycurves(kk).n;
    delta_p = p_cond - p_ctrl;
    err = sqrt((psycurves(kk).err .^ 2) + (psycurves(kk).err .^ 2));
    errorbar(normstim_set,delta_p,err,...
        'color','k',...
        'marker',psyopt.plot.marker,...
        'markersize',psyopt.plot.markersize-2.5,...
        'markerfacecolor',contrast_clrs(kk,:),...
        'linewidth',psyopt.plot.linewidth,...
        'capsize',0);
end

%% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end