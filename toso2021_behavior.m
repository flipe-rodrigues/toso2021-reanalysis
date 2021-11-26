%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% choice GLM (complete model)
X = [t1,t2,i1,i2];
Z = (X - nanmean(X)) ./ nanstd(X);
mdl = fitglm(Z(valid_flags,:),choices(valid_flags,:),'linear',...
    'predictorvars',{'T1','T2','I1','I2'},...
    'distribution','binomial',...
    'intercept',false);
betas = mdl.Coefficients.Estimate;
n_betas = numel(betas);
fig = figure(figopt,...
    'name','choice_GLM');
axes(axesopt.default,...
    'xlim',[0,n_betas+1],...
    'xtick',1:n_betas,...
    'xticklabel',{'T_1','T_2','I_1','I_2'});
title('T_{2}>T_{1} ~ Binomial(\phi(\betaX))');
xlabel('X');
ylabel('\beta');

% plot coefficients
p = stem(1:numel(betas),betas,...
    'color','k',...
    'marker','o',...
    'markersize',10,...
    'markerfacecolor','k',...
    'markeredgecolor','w',...
    'linewidth',1.5);
p.BaseLine.LineWidth = p.LineWidth;

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% choice GLM (partial model)
mdl = fitglm([t1(valid_flags),t2(valid_flags)].*[1,1],...
    choices(valid_flags,:),'linear',...
    'predictorvars',{'T1','T2'},...
    'distribution','binomial',...
    'intercept',true);
betas = mdl.Coefficients.Estimate;
beta_0 = betas(1);
beta_t1 = betas(2);
beta_t2 = betas(3);

%% generalization matrix

% transfer function
tfun = @(x) log(x);

% pair specification
t_pairs = X(:,1:2);
t_pairset = unique(t_pairs(valid_flags,:),'rows');
n_pairs = size(t_pairset,1);

% preallocation
p_choice = nan(n_pairs,1);

% iterate through T1-T2 pairs
for ii = 1 : n_pairs
    t_flags = all(t_pairs == t_pairset(ii,:),2);
    trial_flags = ...
        valid_flags & ...
        t_flags;
    if sum(trial_flags) == 0
        continue;
    end
    
    % compute average performance for the current pair
    p_choice(ii) = mean(choices(trial_flags));
end

% nan filtering
t_pairset = t_pairset(~isnan(p_choice),:);
p_choice = p_choice(~isnan(p_choice));

% figure & axes initialization
fig = figure(figopt,...
    'name','sampling_scheme');
axes(axesopt.default,...
    'xlim',tfun([t_set(1),t_set(end)]) + [-1,1] * .1 * range(tfun(t_set)),...
    'ylim',tfun([t_set(1),t_set(end)]) + [-1,1] * .1 * range(tfun(t_set)),...
    'xtick',tfun(t_set),...
    'ytick',tfun(t_set),...
    'xticklabel',num2cell(t_set),...
    'yticklabel',num2cell(t_set),...
    'xticklabelrotation',45,...
    'yticklabelrotation',45,...
    'xscale','linear',...
    'yscale','linear');
xlabel('T_1 (ms)');
ylabel('T_2 (ms)');

% colorbar
cbar = colorbar(...
    'limits',[.1,.9],...
    'box',axesopt.default.box,...
    'linewidth',axesopt.default.linewidth,...
    'tickdirection','out',...
    'ticklength',unique(axesopt.default.ticklength),...
    'fontsize',axesopt.default.fontsize,...
    'ticks',[.1,.5,.9]);
cbar.Label.String = 'P(T_2 > T_1)';
cbar.Label.Rotation = -90;
cbar.Label.VerticalAlignment = 'bottom';

% plot reference lines
plot(xlim,ylim,'--w',...
    'linewidth',1.5);
plot(tfun([1,1]*median(t1(valid_flags))),ylim,'--w',...
    'linewidth',1.5);

% plot decision boundary
x = linspace(min(t_pairs(:,1)),max(t_pairs(:,1)),1e3);
y = -beta_t1 / beta_t2 * x - beta_0 / beta_t2;
plot(tfun(x),tfun(y),'--k',...
    'linewidth',1.5);

% plot performance
scatter(...
    tfun(t_pairset(:,1)),...
    tfun(t_pairset(:,2)),...
    250,p_choice,'s','filled',...
    'markeredgecolor','k',...
    'linewidth',1.5);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% contraction bias visualization

% figure & axes initialization
fig = figure(figopt,...
    'name','contraction_bias');
axes(axesopt.default,...
    'xlim',tfun([t_set(1),t_set(end)]) + [-1,1] * .1 * range(tfun(t_set)),...
    'ylim',tfun([t_set(1),t_set(end)]) + [-1,1] * .1 * range(tfun(t_set)),...
    'xtick',tfun(t_set),...
    'ytick',tfun(t_set),...
    'xticklabel',num2cell(t_set),...
    'yticklabel',num2cell(t_set),...
    'xticklabelrotation',45,...
    'yticklabelrotation',45,...
    'xscale','linear',...
    'yscale','linear');
xlabel('T_1 (ms)');
ylabel('T_2 (ms)');

% colorbar
cbar = colorbar(...
    'limits',[.1,.9],...
    'box',axesopt.default.box,...
    'linewidth',axesopt.default.linewidth,...
    'tickdirection','out',...
    'ticklength',unique(axesopt.default.ticklength),...
    'fontsize',axesopt.default.fontsize,...
    'ticks',[.1,.5,.9]);
cbar.Label.String = 'P(T_2 > T_1)';
cbar.Label.Rotation = -90;
cbar.Label.VerticalAlignment = 'bottom';

% plot hypothesized probability of reporting T2 > T1
pchoice_x = linspace(min(t1),max(t1),1e3);
beta = 6e-3;
pchoice = 1 - 1 ./ (1 + exp(-beta * (-pchoice_x' + pchoice_x)));
imagesc(xlim,ylim,pchoice);

% plot reference lines
plot(xlim,ylim,'--w',...
    'linewidth',1.5);
plot(tfun([1,1]*median(t1(valid_flags))),ylim,'--w',...
    'linewidth',1.5);

% plot performance
scatter(...
    tfun(t_pairset(:,1)),...
    tfun(t_pairset(:,2)),...
    250,p_choice,'s','filled',...
    'markeredgecolor','none',...
    'linewidth',1.5);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% construct psychophysical triple

% stimulus settings
if strcmpi(task,'duration')
    w1 = 0; % beta_t1;
    w2 = 1; % beta_t2;
    norm_c = (abs(w2) + abs(w1));
%     w2 = w2 / norm_c;
%     w1 = w1 / norm_c;
    stimuli = t2 * w2 + t1 * w1;
else
    w1 = 0;
    w2 = 1;
    norm_c = (abs(w2) + abs(w1));
    w2 = w2 / norm_c;
    w1 = w1 / norm_c;
    stimuli = i2 * w2 + i1 * w1;
end
norm_stimuli = (stimuli - min(stimuli)) / range(stimuli);
stim_set = unique(stimuli(valid_flags));
normstim_set = unique(norm_stimuli(valid_flags));
n_stimuli = numel(stim_set);

% contrast settings
contrast_str = 'i2';
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str(1:end-1),'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);

% preallocation
psycurves = struct();

% iterate through contrasts
for kk = 1 : n_contrasts
    contrast_flags = contrasts == contrast_set(kk);
    
    % iterate through stimuli
    for ii = 1 : n_stimuli
        stim_flags = stimuli == stim_set(ii);
        trial_flags = ...
            valid_flags & ...
            contrast_flags & ...
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
    contrast_flags = contrasts == contrast_set(kk);
    
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
    'name',sprintf('psychometric_curves_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
title(sprintf('Delayed %s comparison task',task));
title('Psychometric curves');
if strcmpi(task,'duration')
    xlabel(sprintf('%.2f \\times T_2 + %.2f \\times T_1 (ms)',w2,w1));
else
    xlabel('I_2 - I_1 (mm/s)');
end
ylabel('P(T_2 > T_1)');

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*median(normstim_set),ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through contrasts
for kk = 1 : n_contrasts
    contrast_flags = contrasts == contrast_set(kk);
    
    % plot psychometric curve
    psyopt.plot.datafaceclr = contrast_clrs(kk,:);
    psyopt.plot.overallvisibility = 'off';
    psyopt.plot.normalizemarkersize = true;
    psyopt.plot.plotfit = true;
    p(kk) = plotpsy(psycurves(kk),psycurves(kk).fit,psyopt.plot);
end

% plot psychometric curve
psyopt.plot.fitclr = 'k';
psyopt.plot.plotfit = true;
psyopt.plot.plotdata = false;

% legend

leg_str = arrayfun(@(x,y)sprintf('I_%s = %i (mm/s)',x,y),...
    repmat(contrast_str(2),n_contrasts,1),contrast_set,...
    'uniformoutput',false);
legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
    'position',[0.085,0.66,.27,.2],...
    'box','on');

%% inset with delta P(long)
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.inset.se,...
    'ylim',[-.25,.25]+[-1,1]*.05*.75,...
    'ytick',[-.25,0,.25]);
ylabel('\DeltaP(T_2 > T_1)',...
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
    if kk == i2_mode_idx
        continue;
    end
    p_ctrl = psycurves(i2_mode_idx).y ./ psycurves(i2_mode_idx).n;
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

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end