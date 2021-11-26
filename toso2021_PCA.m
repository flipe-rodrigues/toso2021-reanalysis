%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% contrast settings
contrast_str = 'i2';
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str(1:end-1),'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);

%% construct S2-aligned, I2-split psths
pre_padd = 500;
roi2use = [0,t_set(end-2)];
roi2plot = [-pre_padd,t_set(end)];
roi2use_n_bins = range(roi2use) * psthbin;
roi2plot_n_bins = range(roi2plot) * psthbin;
roi2use_time = linspace(roi2use(1),roi2use(2),roi2use_n_bins);
roi2plot_time = linspace(roi2plot(1),roi2plot(2),roi2plot_n_bins);
roi2use_flags = ...
    roi2plot_time >= roi2use(1) & ...
    roi2plot_time <= roi2use(2);

% preallocation
s2_psths = nan(roi2plot_n_bins,n_neurons,n_contrasts);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        s2_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(s2_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        s2_spike_counts = data.FR(s2_spike_flags,:);
        s2_spike_rates = conv2(...
            1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
        s2_n_trials = size(s2_spike_counts,1);
        
        % T2-aligned spike rates
        s2_alignment_offset = ...
            pre_init_padding + ...
            pre_t1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            inter_t1t2_delay;
        s2_alignment_flags = ...
            valid_time >= s2_alignment_offset + roi2plot(1) & ...
            valid_time < s2_alignment_offset + t2(s2_spike_flags);
        s2_chunk_flags = ...
            valid_time >= s2_alignment_offset + roi2plot(1) & ...
            valid_time < s2_alignment_offset + roi2plot(2);
        s2_spkrates = s2_spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(...
            s2_spkrates(s2_chunk_flags'),[roi2plot_n_bins,s2_n_trials])';
        
        % compute mean spike density function
        s2_psths(:,nn,ii) = nanmean(s2_spkrates,1);
    end
end

%% normalization
mus = nanmean(s2_psths,[1,3]);
sigs = nanstd(s2_psths,0,[1,3]);
s2_zpsths = (s2_psths - mus) ./ nanstd(sigs);

%% PCA

% concatenate psths across conditions
s2_concat_all = nan(roi2use_n_bins*n_contrasts,n_neurons);
s2_concat_extr = nan(roi2use_n_bins*(n_contrasts-1),n_neurons);
s2_concat_mode = nan(roi2use_n_bins,n_neurons);
s2_concat_diff = nan(roi2use_n_bins,n_neurons);
for nn = 1 : n_neurons
    nn_zpsths_all = s2_zpsths(roi2use_flags,nn,:);
    nn_zpsths_extr = s2_zpsths(roi2use_flags,nn,(1:n_contrasts)~=contrast_mode_idx);
    nn_zpsths_mode = s2_zpsths(roi2use_flags,nn,contrast_mode_idx);
    nn_zpsths_diff = ...
        s2_zpsths(roi2use_flags,nn,end) - s2_zpsths(roi2use_flags,nn,1);
    s2_concat_all(:,nn) = nn_zpsths_all(:);
    s2_concat_extr(:,nn) = nn_zpsths_extr(:);
    s2_concat_mode(:,nn) = nn_zpsths_mode(:);
    s2_concat_diff(:,nn) = nn_zpsths_diff(:);
end

% pca
coeff = pca(s2_concat_extr);
lat_pca = nanvar(s2_concat_all * coeff)';
[~,pca_idcs] = sort(lat_pca,'descend');
coeff = coeff(:,pca_idcs);
exp_pca = lat_pca(pca_idcs) / sum(nanvar(s2_concat_all)) * 100;

% preallocation
s2_score = nan(roi2plot_n_bins,n_neurons,n_contrasts);

% iterate through contrasts
for ii = 1 : n_contrasts
    
    % project onto PCs
    s2_score(:,:,ii) = s2_zpsths(:,:,ii) * coeff;
end

%% 3D trajectories in PC space
fig = figure(figopt,...
    'name','pcTrajectories');
set(gca,...
    axesopt.default,...
    'xtick',0,...
    'ytick',0,...
    'ztick',0);
xlabel(sprintf('%s\n%.1f%% variance','PC 1',exp_pca(1)),...
    'horizontalalignment','center');
ylabel(sprintf('%s\n%.1f%% variance','PC 2',exp_pca(2)),...
    'horizontalalignment','center');
zlabel(sprintf('%s\n%.1f%% variance','PC 3',exp_pca(3)),...
    'horizontalalignment','center');

% iterate through contrasts
for ii = 1 : n_contrasts
    
    % plot trajectory
    plot3(s2_score(:,1,ii),...
        s2_score(:,2,ii),...
        s2_score(:,3,ii),...
        'color',contrast_clrs(ii,:),...
        'linestyle','-',...
        'linewidth',1.5);
    
    % plot stimulus onset
    onset_flags = roi2plot_time <= 0 & ...
        [roi2plot_time(2:end),nan] > 0;
    plot3(s2_score(onset_flags,1,ii),...
        s2_score(onset_flags,2,ii),...
        s2_score(onset_flags,3,ii),...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',5,...
        'markerfacecolor','w',...
        'markeredgecolor',contrast_clrs(ii,:));
    
    % iterate through stimuli
    for tt = 1 : n_t
        
        % plot stimulus offset
        offset_flags = roi2plot_time < t_set(tt) & ...
            [roi2plot_time(2:end),nan] >= t_set(tt);
        plot3(s2_score(offset_flags,1,ii),...
            s2_score(offset_flags,2,ii),...
            s2_score(offset_flags,3,ii),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',6,...
            'markerfacecolor',contrast_clrs(ii,:),...
            'markeredgecolor',contrast_clrs(ii,:));
    end
end

% update axis
axis tight;
xlim(xlim + [-1,1] * .1 * range(xlim));
ylim(ylim + [-1,1] * .1 * range(ylim));
zlim(zlim + [-1,1] * .1 * range(zlim));
set(gca,...
    'xtick',xlim,...
    'ytick',ylim,...
    'ztick',zlim,...
    'xticklabel',{},...
    'yticklabel',{},...
    'zticklabel',{},...
    'xcolor','k',...
    'ycolor','k',...
    'zcolor','k');

% update axis
angle = -135;
view(angle,30);

% save figure
if want2save
    svg_file = fullfile(save_path,'panels',...
        [fig.Name,'_',num2str(angle),'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% PC projections

% figure initialization
fig = figure(figopt,...
    'position',[769.8,41.8,766.4,740.8],...
    'name','pcProjections');
n_pcs2plot = 6;
sps = gobjects(n_pcs2plot,1);
for pc = 1 : n_pcs2plot
    sp_idx = pc * 2 - 1 - (pc > n_pcs2plot / 2) * (n_pcs2plot - 1);
    sps(pc) = subplot(n_pcs2plot/2,2,sp_idx);
    xlabel(sps(pc),'Time since T_2 onset (s)');
    ylabel(sps(pc),sprintf('PC %i\n%.1f%% variance',pc,exp_pca(pc)));
end
xxtick = unique([roi2plot';0;t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
set(sps,...
    axesopt.default,...
    'xlim',roi2plot + [-1,1] * .05 * range(roi2plot),...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'plotboxaspectratio',[2,1,1]);

% link axes
linkaxes(sps,'x');

% iterate through pcs
for pc = 1 : n_pcs2plot
    
    % graphical object preallocation
    p = gobjects(n_contrasts,1);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        
        % plot projection
        p(ii) = plot(sps(pc),roi2plot_time,...
            s2_score(:,pc,ii),...
            'color',contrast_clrs(ii,:),...
            'linestyle','-',...
            'linewidth',1.5);
        
        % plot projection onset
        onset_flags = roi2plot_time <= 0 & ...
            [roi2plot_time(2:end),nan] > 0;
        plot(sps(pc),roi2plot_time(onset_flags),...
            s2_score(onset_flags,pc,ii),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',5,...
            'markerfacecolor','w',...
            'markeredgecolor',contrast_clrs(ii,:));
        
        % iterate through stimuli
        for tt = 1 : n_t
            
            % plot projection offset
            offset_flags = roi2plot_time < t_set(tt) & ...
                [roi2plot_time(2:end),nan] >= t_set(tt);
            plot(sps(pc),roi2plot_time(offset_flags),...
                s2_score(offset_flags,pc,ii),...
                'linewidth',1.5,...
                'marker','o',...
                'markersize',6,...
                'markerfacecolor',contrast_clrs(ii,:),...
                'markeredgecolor','none');
        end
    end
    
    % update axes
    set(sps(pc),...
        'ytick',unique([0,ylim(sps(pc))]),...
        'yticklabel',{'','0',''},...
        'ylim',ylim(sps(pc))+[-1,1]*.1*range(ylim(sps(pc))));
    
    % ui stacking
    uistack(p(:),'bottom');
end

% save figure
if want2save
    svg_file = fullfile(save_path,'panels',[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end