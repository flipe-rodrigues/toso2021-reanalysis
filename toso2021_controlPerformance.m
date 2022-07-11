%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% stimulus settings
w_norm = sum(abs(beta_s1) + abs(beta_s2));
w1 = beta_s1 / w_norm;
w2 = beta_s2 / w_norm;
w1 = -1;
w2 = 1;
stimuli = round(s2 * w2 + s1 * w1);
% stimuli = [ones(size(design,1),1),design] * coeffs;
% % stimuli = [ones(size(Z,1),1),Z] * betas;
% stimuli = round(stimuli,1);
% % stimuli = design(:,7:end) * coeffs(8:end,:);
% valid_flags = ~isnan(stimuli);
stim_set = unique(stimuli(valid_flags));
stim2group_flags = abs(diff(stim_set)) <= range(stim_set) * .01;
if any(stim2group_flags)
    for ii = find(stim2group_flags)'
        stimuli(stimuli == stim_set(ii)) = stim_set(ii + 1);
    end
end
stim_set = unique(stimuli(valid_flags));
n_stimuli = numel(stim_set);
stim_lbl = sprintf('%.2f \\times %s + %.2f \\times %s (%s)',...
    w2,s2_lbl,w1,s1_lbl,s_units);

%% construct psychophysical triple

% normalize stimulus range
norm_stimuli = (stimuli - min(stimuli)) / range(stimuli);
normstim_set = unique(norm_stimuli(valid_flags));

% preallocation
psycurves = struct();

% iterate through contrasts
for kk = 1 : n_contrasts
    contrast_flags = contrasts == contrast_set(kk);
    contrast_flags = i1 == i_set(i1_mode_idx) & i2 == i_set(i2_mode_idx);
    % iterate through stimuli
    for ii = 1 : n_stimuli
        stim_flags = stimuli == stim_set(ii);
        trial_flags = ...
            valid_flags & ...
            ...contrast_flags & ...
            stim_flags;
        psycurves(kk).x(ii,1) = normstim_set(ii);
        psycurves(kk).y(ii,1) = sum(choices(trial_flags));
        psycurves(kk).n(ii,1) = sum(trial_flags);
        psycurves(kk).err(ii,1) = ...
            std(choices(trial_flags)) / sqrt(sum(trial_flags));
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

% iterate through contrasts
for kk = 1 : n_contrasts

    % fit psychometric curve
    psycurves(kk).fit = ...
        psignifit([psycurves(kk).x,psycurves(kk).y,psycurves(kk).n],psyopt.fit);
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
    'name',sprintf('psychometric_curves'));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
title('Psychometric curves');
xlabel(stim_lbl);
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*median(normstim_set),ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through contrasts
for kk = 1 %: n_contrasts

    % plot psychometric curve
    psyopt.plot.datafaceclr = [0,0,0];contrast_clrs(kk,:);
    psyopt.plot.overallvisibility = 'off';
    psyopt.plot.normalizemarkersize = true;
    psyopt.plot.plotfit = false;
    p(kk) = plotpsy(psycurves(kk),psycurves(kk).fit,psyopt.plot);
    plot(psycurves(kk).x,...
        psycurves(kk).y./psycurves(kk).n,...
        '--k','linewidth',1.5);
    plot(psycurves(kk).x(1:n_stimuli/2),...
        psycurves(kk).y(1:n_stimuli/2)./psycurves(kk).n(1:n_stimuli/2),...
        '-k','linewidth',1.5);
    plot(psycurves(kk).x(n_stimuli/2+1:end),...
        psycurves(kk).y(n_stimuli/2+1:end)./psycurves(kk).n(n_stimuli/2+1:end),...
        '-k','linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end