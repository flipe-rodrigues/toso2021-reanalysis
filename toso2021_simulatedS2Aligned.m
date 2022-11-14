%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct T1-aligned, Ii-split psths
pre_padd = 500;
roi = [-pre_padd,t_set(end)];
roi_n_bins = range(roi) * psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

% preallocation
psths = nan(roi2plot_n_bins,n_neurons,n_i);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    i2_flags = i2 == i_set(i2_mode_idx);
    spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i2_flags;
    if sum(spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    spike_counts = data.FR(spike_flags,:);
    spike_rates = conv2(...
        1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
    n_trials = size(spike_counts,1);
    
    % S2-aligned spike rates
    alignment_onset = ...
        pre_init_padding + ...
        pre_t1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi;
    alignment_flags = ...
        valid_time >= alignment_onset + roi2plot(1) & ...
        valid_time < alignment_onset + t2(spike_flags);
    chunk_flags = ...
        valid_time >= alignment_onset + roi2plot(1) & ...
        valid_time < alignment_onset + roi2plot(2);
    spkrates = spike_rates;
    spkrates(~alignment_flags') = nan;
    spkrates = reshape(...
        spkrates(chunk_flags'),[roi2plot_n_bins,n_trials])';
    
    % compute mean spike density function
    psths(:,nn,i2_mode_idx) = nanmean(spkrates,1);
end

%% generate I2-modulated spike rates

% modulation settings
modulation = log(i_set) ./ log(i_set(i2_mode_idx));
scaling = modulation * 1;
gain = modulation * 1;

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'generating fake I2 rates');

    % iterate through intensities
    for ii = 1 : n_i
        if ii == i2_mode_idx
            continue;
        end
        
        % apply temporal scaling
        psths(:,nn,ii) = ...
            interp1(roi_time,psths(:,nn,i2_mode_idx),roi_time*scaling(ii));
    end
end

%% generate I2-modulated spike trains

% preallocation
data.FakeFR = nan(size(data.FR));

% iterate through neurons
for nn = 215 % : n_neurons
    progressreport(nn,n_neurons,'generating fake I2 spikes');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % iterate through intensities
    for ii = 1 : n_i
        i2_flags = i2 == i_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i2_flags;
        spike_trials = find(spike_flags);
        n_trials = numel(spike_trials);
        
        % iterate through trials
        for kk = 1 : n_trials
            trial_idx = spike_trials(kk);
            time_flags = padded_time <= t2(trial_idx);
            lambda = psths(time_flags,nn,i2(trial_idx)==i_set);
            dur = t2(trial_idx) + abs(kernel.paddx(1));
            [n,ts] = poissonprocess(lambda,dur);
            offset = ...
                pre_init_padding + ...
                pre_s1_delay(trial_idx) + ...
                t1(trial_idx) + ...
                isi;
            data.FakeFR(spike_rows(kk),:) = ...
                histcounts(ts+offset,padded_time);
        end
    end
end

%% normalization
mus = nanmean(psths,[1,3]);
sigs = nanstd(psths,0,[1,3]);
s2_zpsths = (psths - mus) ./ sigs;

%% PCA

% concatenate psths across conditions
s2_concat_all = nan(roi_n_bins*n_contrasts,n_neurons);
s2_concat_extr = nan(roi_n_bins*(n_contrasts-1),n_neurons);
s2_concat_mode = nan(roi_n_bins,n_neurons);
s2_concat_diff = nan(roi_n_bins,n_neurons);
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

% training settings
pca_design = s2_concat_all;
%   'all'   ->  vanilla PCA
%   'extr'  ->  pseudo-demixed PCA
%   'mode'  ->  robust PCA

% compute observation weights
weights = ones(size(pca_design,1),1);
for ii = 1 : n_contrasts
    i2_flags = contrasts == contrast_set(ii);
    time_mat = repmat(roi(1) + psthbin : psthbin : roi(2),...
        sum(i2_flags),1);
    concat_idcs = (1 : roi_n_bins) + roi_n_bins * (ii - 1);
    weights(concat_idcs) = sum(time_mat <= t2(i2_flags));
end
weights = repmat(weights,1,size(pca_design,1)/numel(weights));
weights(weights == 0) = nan;

% PCA
coeff = pca(pca_design,...
    'weights',weights);
% coeff_choice = coeff;
% coeff = coeff_choice;

% reorder PCs by variance explained
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
    'name',sprintf('pc_trajectories_t2_%s',contrast_str));
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
        
        % plot cross-condition offset-connecting line
        if ii == n_contrasts
            p = plot3(squeeze(s2_score(offset_flags,1,:)),...
                squeeze(s2_score(offset_flags,2,:)),...
                squeeze(s2_score(offset_flags,3,:)),...
                'linewidth',1,...
                'linestyle','--',...
                'color',contrast_clrs(contrast_mode_idx,:));
            uistack(p,'bottom');
        end
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
angle = 0;
view(angle,0);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% PC projections

% figure initialization
fig = figure(figopt,...
    'position',[535,130,965,860],...
    'name',sprintf('pc_projections_t2_%s',contrast_str));
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
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end