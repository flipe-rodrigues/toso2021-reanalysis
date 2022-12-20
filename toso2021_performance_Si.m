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

%% compute NTD level associated with each stimulus
stim_ntd = nan(n_stimuli,1);

% iterate through stimuli
for ii = 1 : n_stimuli
    if ii <= floor(n_stimuli/2)
        ntd_idx = 1;
    else
        ntd_idx = 2;
    end
    stim_ntd(ii) = ntd_set(ntd_idx);
end

%% NTD colors
ntd_clrs = flipud(gray(n_ntd));

%% construct psychophysical triples

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
            unique_flags & ...
            subject_flags & ...
            stim_flags;
        
        % subject's psychophysical triple
        psycurves(ss).x(ii,1) = normstim_set(ii);
        psycurves(ss).y(ii,1) = sum(choice(trial_flags));
        psycurves(ss).n(ii,1) = sum(trial_flags);
        psycurves(ss).err(ii,1) = ...
            std(choice(trial_flags)) / sqrt(sum(trial_flags));
    end
end

% pooled psychophysical triple
bigpsy.x = normstim_set;
bigpsy.y = sum(horzcat(psycurves.y),2);
bigpsy.n = sum(horzcat(psycurves.n),2);
bigpsy.err = zeros(n_stimuli,1);

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
    
    % fit subject's psychometric curve
    psycurves(ss).fit = ...
        psignifit([psycurves(ss).x,psycurves(ss).y,psycurves(ss).n],psyopt.fit);
end

% fit pooled psychometric curve
bigpsy.fit = psignifit([bigpsy.x,bigpsy.y,bigpsy.n],psyopt.fit);

%% compute stimulus- & NTD-split performance

% preallocation
stim_subj_perf = nan(n_stimuli,n_subjects);
stim_perf = nan(n_stimuli,1);
ntd_subj_perf = nan(n_ntd,n_subjects);
ntd_perf = nan(n_ntd,1);

% iterate through NTD levels
for ii = 1 : n_ntd
    ntd_flags = ...
        valid_flags & ...
        unique_flags & ...
        ntd == ntd_set(ii);
    
    % iterate through subjects
    for jj = 1 : n_subjects
        subject_flags = ...
            ntd_flags & ...
            subjects == subject_set(jj);
        
        % compute performance
        ntd_subj_perf(ii,jj) = ...
            sum(choice(subject_flags)) / sum(subject_flags);
    end
    
    % compute performance
    ntd_perf(ii) = ...
        sum(choice(ntd_flags)) / sum(ntd_flags);
end


% iterate through stimuli
for ii = 1 : n_stimuli
    stim_flags = ...
        valid_flags & ...
        stimuli == stim_set(ii);
    
    % iterate through subjects
    for jj = 1 : n_subjects
        subject_flags = ...
            stim_flags & ...
            subjects == subject_set(jj);
        
        % compute performance
        stim_subj_perf(ii,jj) = ...
            nansum(choice(subject_flags)) / nansum(subject_flags);
    end
    
    % compute performance
    stim_perf(ii) = ...
        sum(choice(stim_flags)) / sum(stim_flags);
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

% figure initialization
fig = figure(figopt,...
    'name',sprintf('performance_Si'));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
xlabel(stim_lbl);
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% reference lines
plot([1,1]*median(normstim_set),ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through subjects
for ss = 1 : n_subjects
    offset = (ss - (n_subjects + 1) / 2) * .0 * range(xlim);
    
    % NTD-crossing line
    plot(normstim_set+offset,stim_subj_perf(:,ss),...
        'color',subject_clr,...
        'linestyle','--',...
        'linewidth',1.5);
    
    % iterate through NTD levels
    for kk = 1 : n_ntd
        ntd_flags = stim_ntd == ntd_set(kk);
        
        % plot subject's performance
        plot(normstim_set(ntd_flags)+offset,...
            stim_subj_perf(ntd_flags,ss),...
            'color',subject_clr,...
            'marker','o',...
            'markersize',6.5,...
            'markeredgecolor',subject_clr,...
            'markerfacecolor',subject_clr,...
            'linestyle','-',...
            'linewidth',1.5);
    end
end

% NTD-crossing line
plot(normstim_set,stim_perf,...
    'color','k',...
    'linestyle','--',...
    'linewidth',1.5);

% iterate through NTD levels
for kk = 1 : n_ntd
    ntd_flags = stim_ntd == ntd_set(kk);
    
    % plot pooled performance
    plot(normstim_set(ntd_flags),...
        stim_perf(ntd_flags),...
        'color','k',...
        'marker','o',...
        'markersize',8.5,...
        'markeredgecolor','k',...
        'markerfacecolor',ntd_clrs(kk,:),...
        'linestyle','-',...
        'linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot performance as a function of stimulus & NTD

% figure initialization
fig = figure(figopt,...
    'position',[744 630 560 412.5],...
    'name',sprintf('performance_NTD_distros'));

% axes initialization
axes(...
    axesopt.default,...
    'plotboxaspectratio',[1,2.25,1],...
    'xlim',[min(ntd_set),max(ntd_set)]+[-1,1]*1/2*range(ntd_set),...
    'xtick',ntd_set,...
    'ylim',axesopt.psycurve.ylim,...
    'ytick',axesopt.psycurve.ytick,...
    'yticklabel',axesopt.psycurve.yticklabel);
xlabel('NTD');
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% reference lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through subjects
for ss = 1 : n_subjects
    offset = (ss - (n_subjects + 1) / 2) * .025 * range(xlim);
    
    % NTD-crossing line
    plot(ntd_set+offset,...
        psycurves(ss).y(floor(n_stimuli/2)+(0:n_ntd-1))./...
        psycurves(ss).n(floor(n_stimuli/2)+(0:n_ntd-1)),...
        'color',subject_clr,...
        'linestyle','--',...
        'linewidth',1.5);
    
    % iterate through NTD levels
    for kk = 1 : n_ntd
        ntd_flags = stim_ntd == ntd_set(kk);
        
        % plot subject's performance
        plot(stim_ntd(ntd_flags)+offset,...
            stim_subj_perf(ntd_flags,ss),...
            'color',subject_clr,...
            'marker','o',...
            'markersize',6.5,...
            'markeredgecolor',subject_clr,...
            'markerfacecolor',subject_clr,...
            'linestyle','-',...
            'linewidth',1.5);
    end
end

% NTD-crossing line
plot(ntd_set,...
    bigpsy.y(floor(n_stimuli/2)+(0:n_ntd-1))./...
    bigpsy.n(floor(n_stimuli/2)+(0:n_ntd-1)),...
    'color','k',...
    'linestyle','--',...
    'linewidth',1.5);

% iterate through NTD levels
for kk = 1 : n_ntd
    ntd_flags = stim_ntd == ntd_set(kk);
    
    % plot average performance
    plot(stim_ntd(ntd_flags),...
        bigpsy.y(ntd_flags)./bigpsy.n(ntd_flags),...
        'color','k',...
        'marker','o',...
        'markersize',8.5,...
        'markeredgecolor','k',...
        'markerfacecolor',ntd_clrs(kk,:),...
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
    'position',[744 630 560 412.5],...
    'name',sprintf('performance_NTD_averages'));

% axes initialization
axes(...
    axesopt.default,...
    'plotboxaspectratio',[1,2.25,1],...
    'xlim',[min(ntd_set),max(ntd_set)]+[-1,1]*1/2*range(ntd_set),...
    'xtick',ntd_set,...
    'ylim',axesopt.psycurve.ylim,...
    'ytick',axesopt.psycurve.ytick,...
    'yticklabel',axesopt.psycurve.yticklabel);
xlabel('NTD');
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% reference lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through subjects
for ss = 1 : n_subjects
    offset = (ss - (n_subjects + 1) / 2) * .025 * range(xlim);
    
    % plot subject's average performance
    plot(ntd_set+offset,ntd_subj_perf(:,ss),...
        'color',subject_clr,...
        'marker','s',...
        'markersize',7.65,...
        'markeredgecolor',subject_clr,...
        'markerfacecolor',subject_clr,...
        'linestyle','-',...
        'linewidth',1.5);
end

% plot pooled performance
plot(ntd_set,ntd_perf,...
    'color','k',...
    'markeredgecolor','k',...
    'markerfacecolor','k',...
    'linestyle','-',...
    'linewidth',1.5);

% iterate through NTD levels
for kk = 1 : n_ntd
    
    % plot pooled performance
    plot(ntd_set(kk),ntd_perf(kk),...
        'color','k',...
        'marker','s',...
        'markersize',12,...
        'markeredgecolor','k',...
        'markerfacecolor',ntd_clrs(kk,:),...
        'linestyle','-',...
        'linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end