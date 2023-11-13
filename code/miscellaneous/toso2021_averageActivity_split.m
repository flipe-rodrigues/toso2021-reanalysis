%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% ROI settings

% roi definition
roi = [-500,t_set(end)];
padded_roi = roi + [-1,1] * .05 * range(roi);
n_bins = range(padded_roi) / psthbin;
time = linspace(padded_roi(1),padded_roi(2),n_bins);

% preallocation
ref_psths = nan(n_bins,n_neurons);
psths = nan(n_bins,n_neurons,n_contrasts);
R = nan(n_neurons,n_contrasts);

%% construct s2-aligned psths

% clamping
i1_clamp_flags = i1 == i_set(i1_mode_idx);
i2_clamp_flags = i2 == i_set(i2_mode_idx);
bw = nan(n_neurons,n_contrasts);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    ref_spike_flags = ...
        valid_flags & ...
        neuron_flags;
    if sum(ref_spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    ref_spike_counts = data.FR(ref_spike_flags,:);
    ref_spike_rates = data.SDF(ref_spike_flags,:);
    %     ref_spike_rates = conv2(...
    %         1,kernel.pdf,ref_spike_counts,'valid')' / psthbin * 1e3;
    ref_n_trials = size(ref_spike_counts,1);
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    %     ref_spike_rates = downsamplecounts(ref_spike_counts,min(t_set));
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    % T2-aligned spike rates
    ref_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(ref_spike_flags) + ...
        t1(ref_spike_flags) + ...
        isi;
    ref_alignment_flags = ...
        padded_time >= ref_alignment + padded_roi(1) & ...
        padded_time < ref_alignment + t2(ref_spike_flags);
    ref_chunk_flags = ...
        padded_time >= ref_alignment + padded_roi(1) & ...
        padded_time < ref_alignment + padded_roi(2);
    ref_spkrates = ref_spike_rates';
    ref_spkrates(~ref_alignment_flags') = nan;
    ref_spkrates = reshape(...
        ref_spkrates(ref_chunk_flags'),[n_bins,ref_n_trials])';
    
    % compute mean spike density function
    ref_psths(:,nn) = nanmean(ref_spkrates,1);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            ...i1_clamp_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = data.SDF(spike_flags,:);
        n_trials = size(spike_counts,1);
        
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        %         spike_rates = downsamplecounts(spike_counts,min(t_set));
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        time_mat = repmat(padded_time,n_trials,1);
        spike_times = time_mat(spike_counts >= 1);
        [~,~,bw(nn,ii)] = ksdensity(sort(spike_times(:)));
        
        % T2-aligned spike rates
        alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi;
        alignment_flags = ...
            padded_time >= alignment_onset + padded_roi(1) & ...
            padded_time < alignment_onset + t2(spike_flags);
        chunk_flags = ...
            padded_time >= alignment_onset + padded_roi(1) & ...
            padded_time < alignment_onset + padded_roi(2);
        spkrates = spike_rates';
        spkrates(~alignment_flags') = nan;
        spkrates = reshape(...
            spkrates(chunk_flags'),[n_bins,n_trials])';
        
        % compute mean spike density function
        psths(:,nn,ii) = nanmean(spkrates,1);
        R(nn,ii) = nanmean(spkrates(:,time>=-334&time<=0),[1,2]);
    end
end

%% normalization

% z-scoring
mus = nanmean(ref_psths,1);
sigs = nanstd(ref_psths,0,1);
zpsths = (psths - mus) ./ sigs;

%% plot overall modulation

% figure initialization
fig = figure(figopt,...
    'name',['average_activity_split_',contrast_str]);

% axes initialization
xxtick = unique([0;roi';t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
sps = gobjects(2,1);
sps(1) = axes(axesopt.default,...
    'xlim',roi + [-1,2] * .05 * range(roi),...
    'xcolor','r',...
    'ycolor','r');
sps(2) = axes(axesopt.default,...
    'xlim',roi + [-2,1] * .05 * range(roi));
set(sps,...
    'clipping','off',...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylim',[-.75,.5]+[-1,1]*.05*1.25,...
    'ytick',[-.75,0,.5]);

xlabel(sps(1),'Time since S_2 onset (ms)');
xlabel(sps(2),'Time since S_2 onset (ms)');
ylabel(sps(1),'Firing rate (z-scored)');
ylabel(sps(2),'Firing rate (z-scored)');

% zero lines
plot([1,1]*0,ylim,...
    'color','k',...
    'linestyle',':');
plot(xlim,[1,1]*0,...
    'color','k',...
    'linestyle',':');

% compute observation weights
time_mat = repmat(...
    padded_roi(1) + psthbin : psthbin : padded_roi(2),n_total_trials,1);
weights = sum(time_mat <= t2);
weights = (weights - min(weights)) ./ range(weights);
t2offset_idcs = find([diff(weights) ~= 0, true]);

% snippet settings
time_snippets = [...
    padded_roi(1),0;...
    0,padded_roi(2)];
n_snippets = size(time_snippets,2);

% iterate through snippets
for ss = 1 : n_snippets
    snippet_flags = ...
        time >= time_snippets(1,ss) & ...
        time <= time_snippets(2,ss);
    
    % graphical object preallocation
    p = gobjects(n_contrasts,1);
    mu_gobjs = gobjects(n_contrasts,n_t);
    sem_gobjs = gobjects(n_contrasts,n_t);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        
        % compute modulation stats
        nan_flags = all(isnan(psths(:,:,ii)),2);
        time_flags = ...
            snippet_flags' & ...
            ~nan_flags;
        X = zpsths(:,:,ii);
        x_mu = nanmean(X,2);
        x_sig = nanstd(X,0,2);
        x_sem = x_sig / sqrt(n_neurons);
        
        % patch s.e.m.
        sem_xpatch = [time(time_flags),fliplr(time(time_flags))];
        sem_ypatch = [x_mu(time_flags)-x_sem(time_flags);...
            flipud(x_mu(time_flags)+x_sem(time_flags))]';
        sem_apatch = [weights(time_flags),fliplr(weights(time_flags))];
        sem_apatch = sem_apatch * range(alphabounds_sem) + alphabounds_sem(1);
        if fadeifnoisy
            alpha_levels = unique(sem_apatch,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = sem_apatch == alpha_levels(aa);
                sem_gobjs(ii,tt) = patch(sps(ss),...
                    sem_xpatch(alpha_flags),...
                    sem_ypatch(alpha_flags),0,...
                    'facealpha',alpha_levels(aa),...
                    'edgecolor','none',...
                    'facecolor',contrast_clrs(ii,:));
            end
        else
            sem_gobjs(ii,1) = patch(sem_xpatch,sem_ypatch,contrast_clrs(ii,:),...
                'facealpha',.25,...
                'edgecolor','none');
        end
        
        % patch average activity
        mu_xpatch = time(time_flags);
        mu_ypatch = x_mu(time_flags)';
        mu_ypatch(end) = nan;
        mu_apatch = weights; % (time_flags);
        mu_apatch = mu_apatch * range(alphabounds_mu) + alphabounds_mu(1);
        if fadeifnoisy
            for tt = 1 : n_t
                try
                    alpha_level = mu_apatch(t2offset_idcs(tt));
                    alpha_flags = mu_apatch(time_flags) == alpha_level;
                    mu_gobjs(ii,tt) = patch(sps(ss),...
                        [mu_xpatch(alpha_flags),nan],...
                        [mu_ypatch(alpha_flags),nan],0,...
                        'edgealpha',alpha_level,...
                        'edgecolor',contrast_clrs(ii,:),...
                        'facecolor','none',...
                        'linewidth',1.5);
                catch
                end
            end
        else
            mu_gobjs(ii,1) = plot(mu_xpatch,mu_ypatch,...
                'color',contrast_clrs(ii,:),...
                'linewidth',1.5);
        end
        
        % plot alignment onset
        %         if ss == 2
        snippet_time = time(time_flags);
        snippet_mu = x_mu(time_flags);
        [~,onset_idx] = min(abs(snippet_time - 0));
        p(ii) = plot(sps(ss),...
            snippet_time(onset_idx),snippet_mu(onset_idx),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',6.5,...
            'markerfacecolor','w',...
            'markeredgecolor',contrast_clrs(ii,:));
        
        % iterate through stimuli
        if ss == 2
            for tt = 1 : n_t
                t2_flags = t2 == t_set(tt);
                trial_flags = ...
                    valid_flags & ...
                    contrast_flags & ...
                    t2_flags;
                if sum(trial_flags) < 1
                    continue;
                end
                
                % plot stimulus offset
                try
                    alpha_level = mu_apatch(t2offset_idcs(tt));
                    scatter(time(t2offset_idcs(tt)),x_mu(t2offset_idcs(tt)),50,...
                        'linewidth',1.5,...
                        'marker','o',...
                        'markerfacealpha',alpha_level.^fadeifnoisy,...
                        'markerfacecolor',contrast_clrs(ii,:),...
                        'markeredgecolor','none');
                catch
                end
            end
        end
    end
    
    % ui restacking
    uistack([mu_gobjs(isgraphics(mu_gobjs)); ...
        sem_gobjs(isgraphics(sem_gobjs))],'bottom');
end

% ui restacking
uistack(p,'top');

% legend
p1 = plot([1,1]*max(xlim)*2,[1,1]*max(ylim)*2,...
    'linewidth',1.5,...
    'linestyle','none',...
    'marker','o',...
    'markersize',7.5,...
    'markerfacecolor','w',...
    'markeredgecolor','k');
p2 = plot([1,1]*max(xlim)*2,[1,1]*max(ylim)*2,...
    'linewidth',1.5,...
    'linestyle','none',...
    'marker','o',...
    'markersize',7.5,...
    'markerfacecolor','k',...
    'markeredgecolor','none');
legend([p1,p2],{'S_2 onset','S_2 offset'},...
    'fontsize',axesopt.default.fontsize*.9,...
    'location','southeast',...
    'box','off');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%%
return;
figure; hold on;
for ii = 1 : n_contrasts
    histogram(bw(:,ii),...
        'facecolor',contrast_clrs(ii,:),...
        'facealpha',.5);
end