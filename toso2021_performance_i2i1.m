%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% stimulus settings
stimuli = round(d2 - d1);
stim_set = unique(stimuli(valid_flags));
stim2group_flags = abs(diff(stim_set)) <= range(stim_set) * .01;
if any(stim2group_flags)
    for ii = find(stim2group_flags)'
        stimuli(stimuli == stim_set(ii)) = stim_set(ii + 1);
    end
end
stim_set = unique(stimuli(valid_flags));
n_stimuli = numel(stim_set);
stim_lbl = sprintf('%s - %s (%s)',d2_lbl,d1_lbl,d_units);

%% compute NID level associated with each stimulus
stim_nid = nan(n_stimuli,1);

% iterate through stimuli
for ii = 1 : n_stimuli
    if ii <= floor(n_stimuli/2)
        nid_idx = 1;
    elseif ii == floor(n_stimuli/2) + 1
        nid_idx = 2;
    else
        nid_idx = 3;
    end
    stim_nid(ii) = nid_set(nid_idx);
end

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

% average psychophysical triple
avgpsy.x = normstim_set;
avgpsy.y = sum(horzcat(psycurves.y),2);
avgpsy.n = sum(horzcat(psycurves.n),2);
avgpsy.err = zeros(n_stimuli,1);

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

% fit average psychometric curve
avgpsy.fit = psignifit([avgpsy.x,avgpsy.y,avgpsy.n],psyopt.fit);

%% compute stimulus- & NID-split performance

% preallocation
nid_stim_perf = cell(n_nid,n_subjects);

% iterate through stimuli
for ii = 1 : n_stimuli
    nid_flags = stim_nid(ii) == nid_set;
    
    % iterate through subjects
    for ss = 1 : n_subjects
        
        % compute performance per stimulus & NID level
        nid_stim_perf{nid_flags,ss} = [nid_stim_perf{nid_flags,ss};...
            psycurves(ss).y(ii) / psycurves(ss).n(ii)];
    end
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
    'name',sprintf('performance_i2i1'));

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

    % plot subject's psychometric curve
    psyopt.plot.datafaceclr = subject_clr;
    psyopt.plot.dataedgeclr = subject_clr;
    psyopt.plot.markersize = 6.5;
    psyopt.plot.overallvisibility = 'off';
    psyopt.plot.normalizemarkersize = false;
    psyopt.plot.plotfit = false;
    plotpsy(psycurves(ss),psycurves(ss).fit,psyopt.plot);
    plot(psycurves(ss).x,...
        psycurves(ss).y./psycurves(ss).n,...
        'color',subject_clr,...
        'linestyle','--',...
        'linewidth',1.5);
    plot(psycurves(ss).x(1:n_stimuli/2),...
        psycurves(ss).y(1:n_stimuli/2)./psycurves(ss).n(1:n_stimuli/2),...
        'color',subject_clr,...
        'linestyle','-',...
        'linewidth',1.5);
    plot(psycurves(ss).x(n_stimuli/2+1:end),...
        psycurves(ss).y(n_stimuli/2+1:end)./psycurves(ss).n(n_stimuli/2+1:end),...
        'color',subject_clr,...
        'linestyle','-',...
        'linewidth',1.5);
end

% plot average psychometric curve
psyopt.plot.datafaceclr = 'k';
psyopt.plot.dataedgeclr = 'k';
psyopt.plot.markersize = 8.5;
psyopt.plot.overallvisibility = 'off';
psyopt.plot.normalizemarkersize = false;
psyopt.plot.plotfit = false;
plotpsy(avgpsy,avgpsy.fit,psyopt.plot);
plot(avgpsy.x,...
    avgpsy.y./avgpsy.n,...
    'color','k',...
    'linestyle','--',...
    'linewidth',1.5);
plot(avgpsy.x(1:n_stimuli/2),...
    avgpsy.y(1:n_stimuli/2)./avgpsy.n(1:n_stimuli/2),...
    'color','k',...
    'linestyle','-',...
    'linewidth',1.5);
plot(avgpsy.x(n_stimuli/2+1:end),...
    avgpsy.y(n_stimuli/2+1:end)./avgpsy.n(n_stimuli/2+1:end),...
    'color','k',...
    'linestyle','-',...
    'linewidth',1.5);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot performance as a function of stimulus & NID

% figure initialization
fig = figure(figopt,...
    'name',sprintf('performance_NID_distros'));

% axes initialization
axes(...
    axesopt.default,...
    'plotboxaspectratio',[1,2.25,1],...
    'xlim',[min(nid_set),max(nid_set)]+[-1,1]*1/2*range(nid_set),...
    'xtick',nid_set,...
    'ylim',axesopt.psycurve.ylim,...
    'ytick',axesopt.psycurve.ytick,...
    'yticklabel',axesopt.psycurve.yticklabel);
xlabel('NID');
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% reference lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through NID levels
for kk = 1 : n_nid
    nid_flags = stim_nid == nid_set(kk);
    
    % iterate through subjects
    for ss = n_subjects : -1 : 1
        offset = (ss - (n_subjects + 1) / 2) * .025 * range(xlim);
        
        % plot subject's performance
        plot(stim_nid(nid_flags)+offset,...
            nid_stim_perf{kk,ss},...
            'color',subject_clr,...
            'marker','o',...
            'markersize',6.5,...
            'markeredgecolor',subject_clr,...
            'markerfacecolor',subject_clr,...
            'linestyle','-',...
            'linewidth',1.5);
        plot(nid_set+offset,...
            psycurves(ss).y(floor(n_stimuli/2)+(0:n_nid-1))./...
            psycurves(ss).n(floor(n_stimuli/2)+(0:n_nid-1)),...
            'color',subject_clr,...
            'linestyle','--',...
            'linewidth',1.5);
    end
end

% iterate through NID levels
for kk = 1 : n_nid
    nid_flags = stim_nid == nid_set(kk);
    
    % plot average performance
    plot(stim_nid(nid_flags)+offset,...
        mean([nid_stim_perf{kk,:}],2),...
        'color','k',...
        'marker','o',...
        'markersize',8.5,...
        'markeredgecolor','k',...
        'markerfacecolor','k',...
        'linestyle','-',...
        'linewidth',1.5);
    plot(nid_set,...
        avgpsy.y(floor(n_stimuli/2)+(0:n_nid-1))./...
        avgpsy.n(floor(n_stimuli/2)+(0:n_nid-1)),...
        'color','k',...
        'linestyle','--',...
        'linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot performance as a function of NID

% figure initialization
fig = figure(figopt,...
    'name',sprintf('performance_NID_averages'));

% axes initialization
axes(...
    axesopt.default,...
    'plotboxaspectratio',[1,2.25,1],...
    'xlim',[min(nid_set),max(nid_set)]+[-1,1]*1/2*range(nid_set),...
    'xtick',nid_set,...
    'ylim',axesopt.psycurve.ylim,...
    'ytick',axesopt.psycurve.ytick,...
    'yticklabel',axesopt.psycurve.yticklabel);
xlabel('NID');
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% reference lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through subjects
for ss = 1 : n_subjects
    offset = (ss - (n_subjects + 1) / 2) * .025 * range(xlim);
    
    % compute average subject's performance per NID level
    nid_mus = cellfun(@(x)mean(x),nid_stim_perf(:,ss));
    nid_std = cellfun(@(x)std(x),nid_stim_perf(:,ss));
    
    % plot subject's average performance
    errorbar(nid_set+offset,nid_mus,nid_std,...
        'color',subject_clr,...
        'marker','s',...
        'markersize',7.65,...
        'markeredgecolor',subject_clr,...
        'markerfacecolor',subject_clr,...
        'linestyle','-',...
        'linewidth',1.5,...
        'capsize',0);
end

% compute average performance per NID level
nid_mus = mean(cellfun(@(x)mean(x),nid_stim_perf),2);
nid_std = std(cellfun(@(x)mean(x),nid_stim_perf),0,2);
nid_sem = nid_std ./ sqrt(n_subjects);

% plot average performance
errorbar(nid_set,nid_mus,nid_sem,...
    'color','k',...
    'marker','s',...
    'markersize',10,...
    'markeredgecolor','k',...
    'markerfacecolor','k',...
    'linestyle','-',...
    'linewidth',1.5,...
    'capsize',0);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end