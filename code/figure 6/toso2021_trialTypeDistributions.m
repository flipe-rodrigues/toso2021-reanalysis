%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% construct T2- & contrast-split trial distributions

% preallocation
trial_counts = nan(n_neurons_total,n_contrasts,n_t);
surviving_trial_counts = nan(n_neurons,n_contrasts,n_t);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);

    % iterate through contrasts
    for kk = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(kk);

        % iterate through S2 durations
        for tt = 1 : n_t
            t2_flags = t2 == t_set(tt);
            trial_flags = ...
                valid_flags & ...
                neuron_flags & ...
                contrast_flags & ...
                t2_flags;

            % compute trial counts
            trial_counts(nn,kk,tt) = sum(trial_flags,'omitnan');
        end

         % compute surviving trial counts
        surviving_trial_counts(nn,kk,:) = ...
            cumsum(trial_counts(nn,kk,:),'reverse','omitnan');
    end
end

%% plot trial type distributions

% figure & axes initialization
fig = figure(figopt,...
    'name',sprintf('trial_type_distributions_%s',contrast_str),...
    'position',[440,230,640,500]);
axes(axesopt.default,...
    'xlim',[0,60],...
    'xtick',[1,2,3,5,10,20,30,50],...
    'ylim',[1,n_t] + [-1,1] * .05 * n_t,...
    'zlim',[0,150],...
    'ztick',linspace(0,150,4),...
    'ytick',1:n_t,...
    'yticklabel',num2cell(t_set),...
    'xscale','log');
xlabel('Trial count',...
    'rotation',-30);
ylabel('T2 (ms)',...
    'rotation',60);
zlabel('Neuron count');
view(30,80);

% bin settings
n_bins = 50;
x_edges = linspace(0,50,n_bins+1);
y_edges = linspace(0,n_t+1,n_bins+1);

% iterate through contrasts
for kk = 1 : n_contrasts

    % iterate through S2 durations
    for tt = 1 : n_t
        
        % compute trial counts
        bin_counts = histcounts2(squeeze(trial_counts(:,kk,tt)),...
            repmat(tt+kk/10,size(trial_counts,1),1),...
            'xbinedges',x_edges,...
            'ybinedges',y_edges);
        
        % plot trial counts
        histogram2('xbinedges',x_edges,...
            'ybinedges',y_edges,...
            'bincounts',bin_counts,...
            'edgecolor','none',...
            'facecolor',contrast_clrs(kk,:),...
            'facealpha',1);
    end
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot surviving trial counts (onset- & offset-aligned)
ti_padd = [-0,0];

% figure initialization
fig = figure(figopt,...
    'position',[165,135,500,200],...
    'name',sprintf('surviving_trial_count_onoff_%s',contrast_str));

% axes initialization
sps = gobjects(2,1);
sps(1) = subplot(1,2,1);
sps(2) = subplot(1,2,2);
set(sps,...
    axesopt.default,...
    'plotboxaspectratiomode','auto',...
    'plotboxaspectratio',[1.1,1,1],...
    'ticklength',axesopt.default.ticklength*2*1.2,...
    'clipping','off',...
    'ylim',[1,200],...
    'ytick',[1,2,5,10,20,50,100,200],...
    'yscale','log');
xxtick = unique([ti_padd(1);-pre_s1_delay;0;t_set;t_set(end)+ti_padd(2)]);
xxtick(xxtick == 0) = abs(xxtick(xxtick == 0));
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,[0,1e3,t_set(t1_mode_idx)])) = {''};
set(sps(1),...
    'xlim',[0,max(t_set)]+ti_padd,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel);
xxtick = unique([-ti_padd(1);0;-t_set;-(t_set(end)+ti_padd(2))]);
xxtick(xxtick == 0) = abs(xxtick(xxtick == 0));
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,-[0,1e3,t_set(t1_mode_idx)])) = {''};
set(sps(2),...
    'xlim',sort(-([0,max(t_set)]+ti_padd)),...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'yaxislocation','right');
xlabel(sps(1),'Time since S_2 onset (ms)');
xlabel(sps(2),'Time since S_2 offset (ms)');
ylabel(sps(1),'Trial count');
ylabel(sps(2),'Trial count',...
    'color','none');

% graphical object preallocation
p_on = gobjects(n_contrasts,1);
p_off = gobjects(n_contrasts,1);

% iterate through contrasts
for kk = 1 : n_contrasts

    % plot surviving trial counts
    counts = squeeze(surviving_trial_counts(:,kk,:));
    counts = [counts,counts(:,end)];
    avg = nanmedian(counts);
    sig = std(counts);
    sem = sig ./ sqrt(n_neurons_total);
    iqr = quantile(counts,[.25,.75]) - nanmedian(counts);
    err = iqr;
    
    % stair-patch hack
    xx = [ti_padd(1);sort(repmat(t_set,2,1))];
    xx = xx(1:end-1);
    yy1 = upsample(avg+err(1,:),2)';
    yy1(2:2:end) = yy1(1:2:end-1);
    yy1 = yy1(1:end-2);
    yy2 = upsample(avg+err(2,:),2)';
    yy2(2:2:end) = yy2(1:2:end-1);
    yy2 = yy2(1:end-2);
    xpatch = [xx;flipud(xx)];
    ypatch = [yy1;flipud(yy2)];
    patch(sps(1),xpatch,ypatch,contrast_clrs(kk,:),...
        'facealpha',1/n_contrasts,...
        'facecolor',contrast_clrs(kk,:),...
        'edgecolor','none');
    p_on(kk) = stairs(sps(1),[ti_padd(1);t_set],avg,...
        'color',contrast_clrs(kk,:),...
        'linewidth',1.5);
    
    patch(sps(2),-xpatch,ypatch,contrast_clrs(kk,:),...
        'facealpha',1/n_contrasts,...
        'facecolor',contrast_clrs(kk,:),...
        'edgecolor','none');
    p_off(kk) = stairs(sps(2),-[ti_padd(1);t_set],avg,...
        'color',contrast_clrs(kk,:),...
        'linewidth',1.5);
end

% ui restacking
uistack(p_on,'top');
uistack(p_off,'top');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end