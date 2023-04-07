%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end
% close all;

%% neuron selection

% selected for being good examples of i2-modulation
if strcmpi(task_str,'duration')
    neurons2plot = [...
        21,24,35,38,62,65,68,72,100,130,205,206,215,224,...
        234,241,356,381,391,393,397,402,406,428,441,448,...
        459,461,462,470,473,493,526,544,553,555,566];
    %     neurons2plot = [...
    %         38,72,205,215,224,391,393,397,402,448,459,462,470,526,566];
    neurons2plot = [...
        393,215,72,459,526];
elseif strcmpi(task_str,'intensity')
    neurons2plot = [...
        19,22,30,61,66,70,100,111,112,115,...
        166,238,243,260,344,408,410];
end
% neurons2plot = flagged_neurons;
% neurons2plot = neuron_idcs;
n_neurons2plot = numel(neurons2plot);

%% version-dependent spike marker
version = ver('matlab');
if contains(version.Release,'2019')
    spike_marker = '.';
    spike_markersize = 2.5;
else
    spike_marker = '.';
    spike_markersize = 2.5;
end

%% construct Si-aligned, Ti- & Ii-split psths
ti_padd = [-500,0];

% clamping
i1_clamp_flags = i1 == i_set(i1_mode_idx);
i2_clamp_flags = i2 == i_set(i2_mode_idx);

% iterate through neurons
for nn = 1 : n_neurons2plot
    progressreport(nn,n_neurons2plot,'parsing neural data');
    neuron_flags = data.NeuronNumb == neurons2plot(nn);
    
    % figure initialization
    fig = figure(figopt,...
        'position',[25,350,1975,375],...
        'name',sprintf('neuron_%i_%i',neurons2plot(nn),fadeifnoisy));
    n_rows = 2;
    n_cols = 4;
    n_sps = (n_rows - 1) * n_cols;
    sps = gobjects(n_sps,1);
    for ii = 1 : n_cols
        sps(ii) = subplot(n_rows,n_cols,ii);
        sps(ii+n_cols) = subplot(n_rows,n_cols,n_cols+ii);
    end
    xxlim = [0,max(t_set)]+ti_padd;
    xxtick = unique([ti_padd(1);-pre_s1_delay;0;t_set;t_set(end)+ti_padd(2)]);
    xxticklabel = num2cell(xxtick);
    xxticklabel(xxtick > 0 & xxtick < 1e3) = {''};
    set(sps,...
        axesopt.default,...
        'layer','top',...
        'plotboxaspectratiomode','auto',...
        'ylimspec','tight');
    set(sps(1+[0,n_cols]),...
        'xlim',xxlim,... + [-1,1] * .05 * range(xxlim),...
        'xtick',xxtick,...
        'xticklabel',xxticklabel);
    set(sps(2+[0,n_cols]),...
        'xlim',xxlim,... + [-1,1] * .05 * range(xxlim),...
        'xtick',xxtick,...
        'xticklabel',xxticklabel);
    set(sps(3+[0,n_cols]),...
        'xlim',xxlim,... + [-1,1] * .05 * range(xxlim),...
        'xtick',xxtick,...
        'xticklabel',xxticklabel);
    set(sps(4+[0,n_cols]),...
        'xlim',xxlim,... + [-1,1] * .05 * range(xxlim),...
        'xtick',xxtick,...
        'xticklabel',xxticklabel);
    set(sps(1:n_cols),...
        'xcolor','none');
    xlabel(sps(1+n_cols),'Time since S_2 onset (s)');
    xlabel(sps(2+n_cols),'Time since S_2 onset (s)');
    xlabel(sps(3+n_cols),'Time since S_2 onset (s)');
    xlabel(sps(4+n_cols),'Time since S_2 onset (s)');
    ylbl_1 = ylabel(sps(1),'Firing rate (Hz)');
    ylbl_2 = ylabel(sps(2),'Firing rate (Hz)');
    ylbl_3 = ylabel(sps(3),'Firing rate (Hz)');
    ylbl_4 = ylabel(sps(4),'Firing rate (Hz)');
    ylbl_5 = ylabel(sps(1+n_cols),'Trial #');
    ylbl_6 = ylabel(sps(2+n_cols),'Trial #');
    ylbl_7 = ylabel(sps(3+n_cols),'Trial #');
    ylbl_8 = ylabel(sps(4+n_cols),'Trial #');
    
    % preallocation
    i1_n_trial_counter = 0;
    t1_n_trial_counter = 0;
    i2_n_trial_counter = 0;
    choice_n_trial_counter = 0;
    i1_boundaries = nan(n_i,1);
    t1_boundaries = nan(n_t,1);
    i2_boundaries = nan(n_i,1);
    choice_boundaries = nan(n_i,1);
    
    % graphical object preallocation
    t1_lines = gobjects(n_t,n_t);
    i1_lines = gobjects(n_i,n_t);
    i2_lines = gobjects(n_i,n_t);
    choice_lines = gobjects(n_choices,n_t);
    t1_patches = gobjects(n_t,n_t);
    i1_patches = gobjects(n_i,n_t);
    i2_patches = gobjects(n_i,n_t);
    choice_patches = gobjects(n_choices,n_t);
    
    % iterate through durations
    for tt = 1 : n_t
        t1_flags = t1 == t_set(tt);
        t1_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            t1_flags;
        if sum(t1_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        t1_spike_counts = data.FR(t1_spike_flags,:);
        t1_spike_rates = conv2(...
            1,kernel.pdf,t1_spike_counts,'valid')' / psthbin * 1e3;
        t1_n_trials = size(t1_spike_counts,1);
        t1_n_tbins = range(xlim(sps(1))) / psthbin;
        
        % T2-aligned, T1-split spike rates
        t1_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(t1_spike_flags) + ...
            t1(t1_spike_flags) + ...
            isi;
        t1_alignment_flags = ...
            valid_time >= t1_alignment_onset + min(xlim(sps(1))) & ...
            valid_time < t1_alignment_onset + t2(t1_spike_flags);
        t1_chunk_flags = ...
            valid_time >= t1_alignment_onset + min(xlim(sps(1))) & ...
            valid_time < t1_alignment_onset + max(xlim(sps(1)));
        t1_spkrates = t1_spike_rates;
        t1_spkrates(~t1_alignment_flags') = nan;
        t1_spkrates = reshape(...
            t1_spkrates(t1_chunk_flags'),[t1_n_tbins,t1_n_trials])';

        % time selection
        time2plot = min(xlim(sps(1))) + psthbin : psthbin : max(xlim(sps(1)));
        time_flags = time2plot <= max(xlim(sps(1)));
        
        % compute mean spike density function
        t1_mu = nanmean(t1_spkrates(:,time_flags),1);
        t1_std = nanstd(t1_spkrates(:,time_flags),0,1);
        t1_trials_throughtime = sum(~isnan(t1_spkrates),1);
        t1_sem = t1_std ./ sqrt(t1_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(t1_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        t1_xpatch = [flagged_time,fliplr(flagged_time)];
        t1_ypatch = [t1_mu(~nan_flags)-t1_sem(~nan_flags),...
            fliplr(t1_mu(~nan_flags)+t1_sem(~nan_flags))];
        patch_alpha = [t1_trials_throughtime(~nan_flags),...
            fliplr(t1_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(t1_trials_throughtime)) / ...
            range(t1_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_sem) + alphabounds_sem(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                t1_patches(tt,aa) = patch(sps(1),...
                    t1_xpatch(alpha_flags),...
                    t1_ypatch(alpha_flags),0,...
                    'facealpha',alpha_levels(aa),...
                    'edgecolor','none',...
                    'facecolor',t1_clrs(tt,:));
            end
        else
            t1_patches(tt,1) = patch(sps(1),t1_xpatch,t1_ypatch,t1_clrs(tt,:),...
                'facealpha',.3,...
                'edgecolor','none');
        end
        
        % patch average activity
        t1_xpatch = flagged_time(~nan_flags);
        t1_ypatch = t1_mu(~nan_flags);
        t1_ypatch(end) = nan;
        patch_alpha = t1_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(t1_trials_throughtime)) / ...
            range(t1_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_mu) + alphabounds_mu(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                t1_lines(tt,aa) = patch(sps(1),...
                    [t1_xpatch(alpha_flags),nan],...
                    [t1_ypatch(alpha_flags),nan],0,...
                    'edgealpha',alpha_levels(aa),...
                    'edgecolor',t1_clrs(tt,:),...
                    'facecolor','none',...
                    'linewidth',1.5);
            end
        else
            t1_lines(tt,1) = plot(sps(1),t1_xpatch,t1_ypatch,...
                'color',t1_clrs(tt,:),...
                'linewidth',1.5);
        end
        
        % plot stimulus onset
        plot(sps(1),time2plot(onset_flags),t1_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',6,...
            'markerfacecolor','w',...
            'markeredgecolor',t1_clrs(tt,:));
        
        % iterate through stimuli
        pseudo_tt = 1;
        for jj = 1 : n_t
            t2_flags = t2 == t_set(jj);
            trial_flags = ...
                t1_spike_flags & ...
                t2_flags;
            if sum(trial_flags) < 1
                continue;
            end
            
            % plot stimulus offset
            offset_flags = time2plot < t_set(jj) & ...
                [time2plot(2:end),nan] >= t_set(jj);
            scatter(sps(1),time2plot(offset_flags),t1_mu(offset_flags),36,...
                'linewidth',1.5,...
                'marker','o',...
                'markerfacealpha',alpha_levels(pseudo_tt).^fadeifnoisy,...
                'markerfacecolor',t1_clrs(tt,:),...
                'markeredgecolor','none');
            pseudo_tt = pseudo_tt + 1;
        end
        
        % plot T1 raster
        t1_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(t1_spike_flags) + ...
            t1(t1_spike_flags) + ...
            isi);
        t1_trial_idcs = (1 : t1_n_trials)' + t1_n_trial_counter;
        t1_trial_idcs_global = data.Trial(t1_spike_flags);
        t1_trial_mat = repmat(t1_trial_idcs,1,n_paddedtimebins);
        t1_spike_trials = t1_trial_mat(t1_spike_counts >= 1);
        t1_spike_times = t1_time_mat(t1_spike_counts >= 1);
        trial_sorter = [t2(t1_spike_flags),choice(t1_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + t1_n_trial_counter;
        t1_sorted_trials = resorted_idcs(t1_spike_trials - t1_n_trial_counter);
        t1_unsorted_trials = t1_trial_idcs(t1_spike_trials - t1_n_trial_counter);
        plot(sps(1+n_cols),t1_spike_times,t1_spike_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(1))) + [0,1,1,0] .* range(xlim(sps(1))) * .1/3;
        ypatch = [.5,.5,t1_n_trials+.5,t1_n_trials+.5] + t1_n_trial_counter;
        patch(sps(1+n_cols),xpatch,ypatch,t1_clrs(tt,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % update trial counters
        t1_n_trial_counter = t1_n_trial_counter + t1_n_trials;
        t1_boundaries(tt) = t1_n_trials;
    end
    
    % iterate through intensities
    for ii = 1 : n_i
        i1_flags = i1 == i_set(ii);
        i1_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i1_flags;
        if sum(i1_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        i1_spike_counts = data.FR(i1_spike_flags,:);
        i1_spike_rates = conv2(...
            1,kernel.pdf,i1_spike_counts,'valid')' / psthbin * 1e3;
        i1_n_trials = size(i1_spike_counts,1);
        i1_n_tbins = range(xlim(sps(2))) / psthbin;
        
        % S1-aligned, I1-split spike rates
        i1_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(i1_spike_flags) + ...
            t1(i1_spike_flags) + ...
            isi;
        i1_alignment_flags = ...
            valid_time >= i1_alignment_onset + min(xlim(sps(2))) & ...
            valid_time < i1_alignment_onset + t2(i1_spike_flags);
        i1_chunk_flags = ...
            valid_time >= i1_alignment_onset + min(xlim(sps(2))) & ...
            valid_time < i1_alignment_onset + max(xlim(sps(2)));
        i1_spkrates = i1_spike_rates;
        i1_spkrates(~i1_alignment_flags') = nan;
        i1_spkrates = reshape(...
            i1_spkrates(i1_chunk_flags'),[i1_n_tbins,i1_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(2))) + psthbin : psthbin : max(xlim(sps(2)));
        time_flags = time2plot <= max(xlim(sps(2)));
        
        % compute mean spike density function
        i1_mu = nanmean(i1_spkrates(:,time_flags),1);
        i1_std = nanstd(i1_spkrates(:,time_flags),0,1);
        i1_trials_throughtime = sum(~isnan(i1_spkrates),1);
        i1_sem = i1_std ./ sqrt(i1_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(i1_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        i1_xpatch = [flagged_time,fliplr(flagged_time)];
        i1_ypatch = [i1_mu(~nan_flags)-i1_sem(~nan_flags),...
            fliplr(i1_mu(~nan_flags)+i1_sem(~nan_flags))];
        patch_alpha = [i1_trials_throughtime(~nan_flags),...
            fliplr(i1_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(i1_trials_throughtime)) / ...
            range(i1_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_sem) + alphabounds_sem(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                i1_patches(ii,aa) = patch(sps(2),...
                    i1_xpatch(alpha_flags),...
                    i1_ypatch(alpha_flags),0,...
                    'facealpha',alpha_levels(aa),...
                    'edgecolor','none',...
                    'facecolor',i1_clrs(ii,:));
            end
        else
            i1_patches(ii,1) = patch(sps(2),i1_xpatch,i1_ypatch,i1_clrs(ii,:),...
                'facealpha',.3,...
                'edgecolor','none');
        end
        
        % patch average activity
        i1_xpatch = flagged_time(~nan_flags);
        i1_ypatch = i1_mu(~nan_flags);
        i1_ypatch(end) = nan;
        patch_alpha = i1_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(i1_trials_throughtime)) / ...
            range(i1_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_mu) + alphabounds_mu(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                i1_lines(ii,aa) = patch(sps(2),...
                    [i1_xpatch(alpha_flags),nan],...
                    [i1_ypatch(alpha_flags),nan],0,...
                    'edgealpha',alpha_levels(aa),...
                    'edgecolor',i1_clrs(ii,:),...
                    'facecolor','none',...
                    'linewidth',1.5);
            end
        else
            i1_lines(ii,1) = plot(sps(2),i1_xpatch,i1_ypatch,...
                'color',i1_clrs(ii,:),...
                'linewidth',1.5);
        end
        
        % plot stimulus onset
        plot(sps(2),time2plot(onset_flags),i1_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',6,...
            'markerfacecolor','w',...
            'markeredgecolor',i1_clrs(ii,:));
        
        % iterate through stimuli
        for jj = 1 : n_t
            t2_flags = t2 == t_set(jj);
            trial_flags = ...
                i1_spike_flags & ...
                t2_flags;
            if sum(trial_flags) < 1
                continue;
            end
            
            % plot stimulus offset
            offset_flags = time2plot < t_set(jj) & ...
                [time2plot(2:end),nan] >= t_set(jj);
            scatter(sps(2),time2plot(offset_flags),i1_mu(offset_flags),36,...
                'linewidth',1.5,...
                'marker','o',...
                'markerfacealpha',alpha_levels(jj).^fadeifnoisy,...
                'markerfacecolor',i1_clrs(ii,:),...
                'markeredgecolor','none');
        end
        
        % plot I1 raster
        i1_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(i1_spike_flags) + ...
            t1(i1_spike_flags) + ...
            isi);
        i1_trial_idcs = (1 : i1_n_trials)' + i1_n_trial_counter;
        i1_trial_mat = repmat(i1_trial_idcs,1,n_paddedtimebins);
        i1_spike_trials = i1_trial_mat(i1_spike_counts >= 1);
        i1_spike_times = i1_time_mat(i1_spike_counts >= 1);
        trial_sorter = [t2(i1_spike_flags),choice(i1_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + i1_n_trial_counter;
        i1_sorted_trials = resorted_idcs(i1_spike_trials - i1_n_trial_counter);
        plot(sps(2+n_cols),i1_spike_times,i1_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(2))) + [0,1,1,0] .* range(xlim(sps(2))) * .1/3;
        ypatch = [.5,.5,i1_n_trials+.5,i1_n_trials+.5] + i1_n_trial_counter;
        patch(sps(2+n_cols),xpatch,ypatch,i1_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % update trial counters
        i1_n_trial_counter = i1_n_trial_counter + i1_n_trials;
        i1_boundaries(ii) = i1_n_trials;
    end
    
    % iterate through intensities
    for ii = 1 : n_i
        i2_flags = i2 == i_set(ii);
        i2_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i2_flags;
        if sum(i2_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        i2_spike_counts = data.FR(i2_spike_flags,:);
        i2_spike_rates = conv2(...
            1,kernel.pdf,i2_spike_counts,'valid')' / psthbin * 1e3;
        i2_n_trials = size(i2_spike_counts,1);
        i2_n_tbins = range(xlim(sps(3))) / psthbin;
        
        % S2-aligned, I2-split spike rates
        i2_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(i2_spike_flags) + ...
            t1(i2_spike_flags) + ...
            isi;
        i2_alignment_flags = ...
            valid_time >= i2_alignment_onset + min(xlim(sps(3))) & ...
            valid_time < i2_alignment_onset + t2(i2_spike_flags);
        i2_chunk_flags = ...
            valid_time >= i2_alignment_onset + min(xlim(sps(3))) & ...
            valid_time < i2_alignment_onset + max(xlim(sps(3)));
        i2_spkrates = i2_spike_rates;
        i2_spkrates(~i2_alignment_flags') = nan;
        i2_spkrates = reshape(...
            i2_spkrates(i2_chunk_flags'),[i2_n_tbins,i2_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(3))) + psthbin : psthbin : max(xlim(sps(3)));
        time_flags = time2plot <= max(xlim(sps(3)));
        
        % compute mean spike density function
        i2_mu = nanmean(i2_spkrates(:,time_flags),1);
        i2_std = nanstd(i2_spkrates(:,time_flags),0,1);
        i2_trials_throughtime = sum(~isnan(i2_spkrates),1);
        i2_sem = i2_std ./ sqrt(i2_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(i2_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        i2_xpatch = [flagged_time,fliplr(flagged_time)];
        i2_ypatch = [i2_mu(~nan_flags)-i2_sem(~nan_flags),...
            fliplr(i2_mu(~nan_flags)+i2_sem(~nan_flags))];
        patch_alpha = [i2_trials_throughtime(~nan_flags),...
            fliplr(i2_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(i2_trials_throughtime)) / ...
            range(i2_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_sem) + alphabounds_sem(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                i2_patches(ii,aa) = patch(sps(3),...
                    i2_xpatch(alpha_flags),...
                    i2_ypatch(alpha_flags),0,...
                    'facealpha',alpha_levels(aa),...
                    'edgecolor','none',...
                    'facecolor',i2_clrs(ii,:));
            end
        else
            i2_patches(ii,1) = patch(sps(3),i2_xpatch,i2_ypatch,i2_clrs(ii,:),...
                'facealpha',.3,...
                'edgecolor','none');
        end
        
        % patch average activity
        i2_xpatch = flagged_time(~nan_flags);
        i2_ypatch = i2_mu(~nan_flags);
        i2_ypatch(end) = nan;
        patch_alpha = i2_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(i2_trials_throughtime)) / ...
            range(i2_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_mu) + alphabounds_mu(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                i2_lines(ii,aa) = patch(sps(3),...
                    [i2_xpatch(alpha_flags),nan],...
                    [i2_ypatch(alpha_flags),nan],0,...
                    'edgealpha',alpha_levels(aa),...
                    'edgecolor',i2_clrs(ii,:),...
                    'facecolor','none',...
                    'linewidth',1.5);
            end
        else
            i2_lines(ii,1) = plot(sps(3),i2_xpatch,i2_ypatch,...
                'color',i2_clrs(ii,:),...
                'linewidth',1.5);
        end
        
        % plot stimulus onset
        plot(sps(3),time2plot(onset_flags),i2_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',6,...
            'markerfacecolor','w',...
            'markeredgecolor',i2_clrs(ii,:));
        
        % iterate through stimuli
        for jj = 1 : n_t
            t2_flags = t2 == t_set(jj);
            trial_flags = ...
                i2_spike_flags & ...
                t2_flags;
            if sum(trial_flags) < 1
                continue;
            end
            
            % plot stimulus offset
            offset_flags = time2plot < t_set(jj) & ...
                [time2plot(2:end),nan] >= t_set(jj);
            scatter(sps(3),time2plot(offset_flags),i2_mu(offset_flags),36,...
                'linewidth',1.5,...
                'marker','o',...
                'markerfacealpha',alpha_levels(jj).^fadeifnoisy,...
                'markerfacecolor',i2_clrs(ii,:),...
                'markeredgecolor','none');
        end
        
        % plot I2 raster
        i2_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(i2_spike_flags) + ...
            t1(i2_spike_flags) + ...
            isi);
        i2_trial_idcs = (1 : i2_n_trials)' + i2_n_trial_counter;
        i2_trial_mat = repmat(i2_trial_idcs,1,n_paddedtimebins);
        i2_spike_trials = i2_trial_mat(i2_spike_counts >= 1);
        i2_spike_times = i2_time_mat(i2_spike_counts >= 1);
        trial_sorter = [t2(i2_spike_flags),choice(i2_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + i2_n_trial_counter;
        i2_sorted_trials = resorted_idcs(i2_spike_trials - i2_n_trial_counter);
        plot(sps(3+n_cols),i2_spike_times,i2_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(3))) + [0,1,1,0] .* range(xlim(sps(3))) * .1/3;
        ypatch = [.5,.5,i2_n_trials+.5,i2_n_trials+.5] + i2_n_trial_counter;
        patch(sps(3+n_cols),xpatch,ypatch,i2_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % update trial counters
        i2_n_trial_counter = i2_n_trial_counter + i2_n_trials;
        i2_boundaries(ii) = i2_n_trials;
    end
    
    % iterate through choices
    for ii = 1 : n_choices
        choice_flags = choice == choice_set(ii);
        %         choice_flags = choice_correct == choice_correct_set(ii);
        choice_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            choice_flags;
        if sum(choice_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        choice_spike_counts = data.FR(choice_spike_flags,:);
        choice_spike_rates = conv2(...
            1,kernel.pdf,choice_spike_counts,'valid')' / psthbin * 1e3;
        choice_n_trials = size(choice_spike_counts,1);
        choice_n_tbins = range(xlim(sps(3))) / psthbin;
        
        % S2-aligned, I2-split spike rates
        choice_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(choice_spike_flags) + ...
            t1(choice_spike_flags) + ...
            isi;
        choice_alignment_flags = ...
            valid_time >= choice_alignment_onset + min(xlim(sps(3))) & ...
            valid_time < choice_alignment_onset + t2(choice_spike_flags);
        choice_chunk_flags = ...
            valid_time >= choice_alignment_onset + min(xlim(sps(3))) & ...
            valid_time < choice_alignment_onset + max(xlim(sps(3)));
        choice_spkrates = choice_spike_rates;
        choice_spkrates(~choice_alignment_flags') = nan;
        choice_spkrates = reshape(...
            choice_spkrates(choice_chunk_flags'),[choice_n_tbins,choice_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(3))) + psthbin : psthbin : max(xlim(sps(3)));
        time_flags = time2plot <= max(xlim(sps(3)));
        
        % compute mean spike density function
        choice_mu = nanmean(choice_spkrates(:,time_flags),1);
        choice_std = nanstd(choice_spkrates(:,time_flags),0,1);
        choice_trials_throughtime = sum(~isnan(choice_spkrates),1);
        choice_sem = choice_std ./ sqrt(choice_trials_throughtime);
        
        %
        t2countthroughtime = ...
            sum(flagged_time <= t2(valid_flags & neuron_flags));
        t2offset_idcs = find([diff(t2countthroughtime) ~= 0, true]);
        
        % flag current stimulus period
        nan_flags = isnan(choice_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        choice_xpatch = [flagged_time,fliplr(flagged_time)];
        choice_ypatch = [choice_mu(~nan_flags)-choice_sem(~nan_flags),...
            fliplr(choice_mu(~nan_flags)+choice_sem(~nan_flags))];
        patch_alpha = [choice_trials_throughtime(~nan_flags),...
            fliplr(choice_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(choice_trials_throughtime)) / ...
            range(choice_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_sem) + alphabounds_sem(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_t
                alpha_level = patch_alpha(t2offset_idcs(aa));
                alpha_flags = patch_alpha == alpha_level;
                choice_patches(ii,aa) = patch(sps(4),...
                    choice_xpatch(alpha_flags),...
                    choice_ypatch(alpha_flags),0,...
                    'facealpha',alpha_level,...
                    'edgecolor','none',...
                    'facecolor',choice_clrs(ii,:));
            end
        else
            choice_patches(ii,1) = patch(sps(4),choice_xpatch,choice_ypatch,choice_clrs(ii,:),...
                'facealpha',.3,...
                'edgecolor','none');
        end
        
        % patch average activity
        choice_xpatch = flagged_time(~nan_flags);
        choice_ypatch = choice_mu(~nan_flags);
        choice_ypatch(end) = nan;
        patch_alpha = choice_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(choice_trials_throughtime)) / ...
            range(choice_trials_throughtime);
        patch_alpha = patch_alpha * range(alphabounds_mu) + alphabounds_mu(1);
        if fadeifnoisy
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                choice_lines(ii,aa) = patch(sps(4),...
                    [choice_xpatch(alpha_flags),nan],...
                    [choice_ypatch(alpha_flags),nan],0,...
                    'edgealpha',alpha_levels(aa),...
                    'edgecolor',choice_clrs(ii,:),...
                    'facecolor','none',...
                    'linewidth',1.5);
            end
        else
            choice_lines(ii,1) = plot(sps(4),choice_xpatch,choice_ypatch,...
                'color',choice_clrs(ii,:),...
                'linewidth',1.5);
        end
        
        % plot stimulus onset
        plot(sps(4),time2plot(onset_flags),choice_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',6,...
            'markerfacecolor','w',...
            'markeredgecolor',choice_clrs(ii,:));
        
        % iterate through stimuli
        for jj = 1 : n_t
            
            % plot stimulus offset
            alpha_level = patch_alpha(t2offset_idcs(jj));
            offset_flags = time2plot < t_set(jj) & ...
                [time2plot(2:end),nan] >= t_set(jj);
            scatter(sps(4),time2plot(offset_flags),choice_mu(offset_flags),36,...
                'linewidth',1.5,...
                'marker','o',...
                'markerfacealpha',alpha_level.^fadeifnoisy,...
                'markerfacecolor',choice_clrs(ii,:),...
                'markeredgecolor','none');
        end
        
        % plot I2 raster
        choice_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(choice_spike_flags) + ...
            t1(choice_spike_flags) + ...
            isi);
        choice_trial_idcs = (1 : choice_n_trials)' + choice_n_trial_counter;
        choice_trial_mat = repmat(choice_trial_idcs,1,n_paddedtimebins);
        choice_spike_trials = choice_trial_mat(choice_spike_counts >= 1);
        choice_spike_times = choice_time_mat(choice_spike_counts >= 1);
        trial_sorter = [t2(choice_spike_flags),choice(choice_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + choice_n_trial_counter;
        choice_sorted_trials = resorted_idcs(choice_spike_trials - choice_n_trial_counter);
        plot(sps(4+n_cols),choice_spike_times,choice_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(4))) + [0,1,1,0] .* range(xlim(sps(4))) * .1/3;
        ypatch = [.5,.5,choice_n_trials+.5,choice_n_trials+.5] + choice_n_trial_counter;
        patch(sps(4+n_cols),xpatch,ypatch,choice_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % update trial counters
        choice_n_trial_counter = choice_n_trial_counter + choice_n_trials;
        choice_boundaries(ii) = choice_n_trials;
    end
    
    % ui restacking
    uistack([t1_lines(isgraphics(t1_lines)); ...
        t1_patches(isgraphics(t1_patches))],'bottom');
    uistack([i1_lines(isgraphics(i1_lines)); ...
        i1_patches(isgraphics(i1_patches))],'bottom');
    uistack([i2_lines(isgraphics(i2_lines)); ...
        i2_patches(isgraphics(i2_patches))],'bottom');
    uistack([choice_lines(isgraphics(choice_lines)); ...
        choice_patches(isgraphics(choice_patches))],'bottom');
    
    %
    if all(isnan(t1_boundaries)) || ...
            all(isnan(i1_boundaries)) || ...
            all(isnan(i2_boundaries)) || ...
            all(isnan(choice_boundaries))
        continue;
    end
    
    % axes updates
    yylim = [1,nansum(t1_boundaries)];
    yytick = unique([1;cumsum(t1_boundaries,'omitnan')]);
    yyticklabel = num2cell(yytick);
    yyticklabel(~ismember(yytick,yylim)) = {''};
    set(sps(1+n_cols),...
        'ylim',yylim+[-1,1]*.5,...
        'ytick',yytick,...
        'yticklabel',yyticklabel);
    yylim = [1,nansum(i1_boundaries)];
    yytick = unique([1;cumsum(i1_boundaries,'omitnan')]);
    yyticklabel = num2cell(yytick);
    yyticklabel(~ismember(yytick,yylim)) = {''};
    set(sps(2+n_cols),...
        'ylim',yylim+[-1,1]*.5,...
        'ytick',yytick,...
        'yticklabel',yyticklabel);
    yylim = [1,nansum(i2_boundaries)];
    yytick = unique([1;cumsum(i2_boundaries,'omitnan')]);
    yyticklabel = num2cell(yytick);
    yyticklabel(~ismember(yytick,yylim)) = {''};
    set(sps(3+n_cols),...
        'ylim',yylim+[-1,1]*.5,...
        'ytick',yytick,...
        'yticklabel',yyticklabel);
    yylim = [1,nansum(choice_boundaries)];
    yytick = unique([1;cumsum(choice_boundaries,'omitnan')]);
    yyticklabel = num2cell(yytick);
    yyticklabel(~ismember(yytick,yylim)) = {''};
    set(sps(4+n_cols),...
        'ylim',yylim+[-1,1]*.5,...
        'ytick',yytick,...
        'yticklabel',yyticklabel);
    
    % update y-axis limits
    yylim = [min([ylim(sps(1)),ylim(sps(2)),ylim(sps(3))]),...
        max([ylim(sps(1)),ylim(sps(2)),ylim(sps(3))])];
    %     yylim = [0,max(ylim(sps(3)))];
    yylim = 5 * [floor(yylim(1)/5),ceil(yylim(2)/5)];
    if ismember(neurons2plot(nn),[215])
        yylim = [15,65];
    elseif ismember(neurons2plot(nn),[462])
        yylim = [0,60];
    elseif ismember(neurons2plot(nn),[402])
        yylim = [0,50];
    elseif ismember(neurons2plot(nn),[381])
        yylim = [0,40];
    elseif ismember(neurons2plot(nn),[68,205,441])
        yylim = [0,15];
    elseif ismember(neurons2plot(nn),[72])
        yylim = [0,15];
    elseif ismember(neurons2plot(nn),[459])
        yylim = [5,50];
    elseif ismember(neurons2plot(nn),[473])
        yylim = [5,25];
    elseif ismember(neurons2plot(nn),[393])
        yylim = [5,30];
    elseif ismember(neurons2plot(nn),[526])
        yylim = [5,40];
    elseif ismember(neurons2plot(nn),[206])
        yylim = [0,80];
    elseif ismember(neurons2plot(nn),[35])
        yylim = [0,25];
    elseif ismember(neurons2plot(nn),[379])
        yylim = [0,5];
    end
    yylim = max(yylim,[0,1]);
    set(sps(1:n_cols),...
        'ylim',yylim,...+[-1,1]*.1*range(yylim),...
        'ytick',yylim,...
        'clipping','off');
    drawnow;
    
    % update y-labels
    ylbl_5.Position(1) = ylbl_1.Position(1);
    ylbl_6.Position(1) = ylbl_2.Position(1);
    ylbl_7.Position(1) = ylbl_3.Position(1);
    ylbl_8.Position(1) = ylbl_4.Position(1);
    
    % shift psth axes down a bit
    for ii = 1 : n_cols
        set(sps(ii),...
            'position',get(sps(ii),'position')-[0,.085,0,0]);
    end
    
    % save figure
    if want2save
        try
            % save settings
            png_file = fullfile(raster_path,[get(fig,'name'),'.png']);
            print(fig,png_file,'-dpng','-r300','-painters');
            svg_file = fullfile(panel_path,[fig.Name,'.svg']);
            print(fig,svg_file,'-dsvg','-painters');
            close(fig);
        catch
            close(fig);
        end
    else
        pause(1);
        close(fig);
    end
end