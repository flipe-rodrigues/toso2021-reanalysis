%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end
close all;

%% neuron selection

% selected for being good examples of i2-modulation
if strcmpi(task_str,'duration')
    neurons2plot = [...
        38,68,72,215,234,381,391,393,402,436,...
        441,448,459,462,470,473,506,526];
    neurons2plot = fliplr([...
        393,473,215,72,459,526]);
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
        ...'windowstate','maximized',...
        'position',[150,100,1750,460],...
        'name',sprintf('neuron_%i',neurons2plot(nn)));
    n_rows = 2;
    n_cols = 4;
    n_sps = (n_rows - 1) * n_cols;
    sps = gobjects(n_sps,1);
    for ii = 1 : n_cols
        sps(ii) = subplot(n_rows,n_cols,ii);
        sps(ii+n_cols) = subplot(n_rows,n_cols,n_cols+ii);
    end
    xxtick_on = unique([ti_padd(1);-pre_s1_delay;0;t_set;t_set(end)+ti_padd(2)]);
    xxticklabel_on = num2cell(xxtick_on);
    xxticklabel_on(xxtick_on > 0 & xxtick_on < 1e3) = {''};
    xxtick_off = unique([-ti_padd(1);0;-t_set;-(t_set(end)+ti_padd(2))]);
    xxticklabel_off = num2cell(xxtick_off);
    xxticklabel_off(xxtick_off < 0 & xxtick_off > -1e3) = {''};
    set(sps,...
        axesopt.default,...
        'layer','top',...
        'plotboxaspectratiomode','auto',...
        'ylimspec','tight');
    set(sps(1+[0,n_cols]),...
        'xlim',[0,max(t_set)]+ti_padd,...
        'xtick',xxtick_on,...
        'xticklabel',xxticklabel_on);
    set(sps(2+[0,n_cols]),...
        'xlim',sort(-([0,max(t_set)]+ti_padd)),...
        'xtick',xxtick_off,...
        'xticklabel',xxticklabel_off);
    set(sps(3+[0,n_cols]),...
        'xlim',[0,max(t_set)]+ti_padd,...
        'xtick',xxtick_on,...
        'xticklabel',xxticklabel_on);
    set(sps(4+[0,n_cols]),...
        'xlim',sort(-([0,max(t_set)]+ti_padd)),...
        'xtick',xxtick_off,...
        'xticklabel',xxticklabel_off);
    set(sps(1:n_cols),...
        'xcolor','none');
    %     set(sps((1:n_cols)+n_cols),...
    %         'plotboxaspectratio',[1,1,1]);
    xlabel(sps(1+n_cols),'Time since S_1 onset (s)');
    xlabel(sps(2+n_cols),'Time since S_1 offset (s)');
    xlabel(sps(3+n_cols),'Time since S_2 onset (s)');
    xlabel(sps(4+n_cols),'Time since S_2 offset (s)');
    ylbl_1 = ylabel(sps(1),'Firing rate (Hz)');
    ylbl_2 = ylabel(sps(2),'Firing rate (Hz)');
    ylbl_3 = ylabel(sps(3),'Firing rate (Hz)');
    ylbl_4 = ylabel(sps(4),'Firing rate (Hz)');
    ylbl_5 = ylabel(sps(1+n_cols),'Trial #');
    ylbl_6 = ylabel(sps(2+n_cols),'Trial #');
    ylbl_7 = ylabel(sps(3+n_cols),'Trial #');
    ylbl_8 = ylabel(sps(4+n_cols),'Trial #');
    
    % waveform axes
    axwav = axes('position',[.915,.5,.075,.15],...
        'nextplot','add',...
        'xlimspec','tight',...
        'ylimspec','tight',...
        'box','off',...
        'xcolor','none',...
        'ycolor','none');
    
    % preallocation
    s1_n_trial_counter = 0;
    s2_n_trial_counter = 0;
    s1_boundaries = nan(n_i,1);
    s2_boundaries = nan(n_i,1);
    
    % iterate through intensities
    for ii = 1 : n_i
        i1_flags = i1 == i_set(ii);
        s1_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i1_flags;
        if sum(s1_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        s1_spike_counts = data.FR(s1_spike_flags,:);
        s1_spike_rates = conv2(...
            1,kernel.pdf,s1_spike_counts,'valid')' / psthbin * 1e3;
        s1_n_trials = size(s1_spike_counts,1);
        s1on_n_tbins = range(xlim(sps(1))) * psthbin;
        s1off_n_tbins = range(xlim(sps(2))) * psthbin;
        
        % S1-onset-aligned, I1-split spike rates
        s1on_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(s1_spike_flags);
        s1on_alignment_flags = ...
            valid_time >= s1on_alignment_onset + min(xlim(sps(1))) & ...
            valid_time < s1on_alignment_onset + t1(s1_spike_flags);
        s1on_chunk_flags = ...
            valid_time >= s1on_alignment_onset + min(xlim(sps(1))) & ...
            valid_time < s1on_alignment_onset + max(xlim(sps(1)));
        s1on_spkrates = s1_spike_rates;
        s1on_spkrates(~s1on_alignment_flags') = nan;
        s1on_spkrates = reshape(...
            s1on_spkrates(s1on_chunk_flags'),[s1on_n_tbins,s1_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(1))) + psthbin : psthbin : max(xlim(sps(1)));
        time_flags = time2plot <= max(xlim(sps(1)));
        
        % compute mean spike density function
        s1on_mu = nanmean(s1on_spkrates(:,time_flags),1);
        s1on_std = nanstd(s1on_spkrates(:,time_flags),0,1);
        s1on_trials_throughtime = sum(~isnan(s1on_spkrates),1);
        s1on_sem = s1on_std ./ sqrt(s1on_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(s1on_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        s1on_xpatch = [flagged_time,fliplr(flagged_time)];
        s1on_ypatch = [s1on_mu(~nan_flags)-s1on_sem(~nan_flags),...
            fliplr(s1on_mu(~nan_flags)+s1on_sem(~nan_flags))];
        patch_alpha = [s1on_trials_throughtime(~nan_flags),...
            fliplr(s1on_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(s1on_trials_throughtime)) / ...
            range(s1on_trials_throughtime);
        patch_alpha = patch_alpha * .25 + .05;
        patch(sps(1),s1on_xpatch,s1on_ypatch,i1_clrs(ii,:),...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'facealpha','interp',...
            'edgecolor','none');
        
        % patch average activity
        s1on_xpatch = flagged_time;
        s1on_ypatch = s1on_mu(~nan_flags);
        s1on_ypatch(end) = nan;
        patch_alpha = s1on_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(s1on_trials_throughtime)) / ...
            range(s1on_trials_throughtime);
        patch_alpha = patch_alpha * .95 + .05;
        patch(sps(1),s1on_xpatch,s1on_ypatch,0,...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'edgealpha','interp',...
            'edgecolor',i1_clrs(ii,:),...
            'linewidth',1.5);
        
        % plot stimulus onset
        plot(sps(1),time2plot(onset_flags),s1on_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',7.5,...
            'markerfacecolor','w',...
            'markeredgecolor',i1_clrs(ii,:));
        
        % plot I1 raster
        s1on_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(s1_spike_flags));
        s1on_trial_idcs = (1 : s1_n_trials)' + s1_n_trial_counter;
        s1on_trial_mat = repmat(s1on_trial_idcs,1,n_paddedtimebins);
        s1on_spike_trials = s1on_trial_mat(s1_spike_counts >= 1);
        s1on_spike_times = s1on_time_mat(s1_spike_counts >= 1);
        trial_sorter = [t1(s1_spike_flags),prev_choice(s1_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + s1_n_trial_counter;
        s1on_sorted_trials = resorted_idcs(s1on_spike_trials - s1_n_trial_counter);
        plot(sps(1+n_cols),s1on_spike_times,s1on_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(1))) + [0,1,1,0] .* range(xlim(sps(1))) * .1/3;
        ypatch = [.5,.5,s1_n_trials+.5,s1_n_trials+.5] + s1_n_trial_counter;
        patch(sps(1+n_cols),xpatch,ypatch,i1_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % S1-offset-aligned, I1-split spike rates
        s1off_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(s1_spike_flags) + ...
            t1(s1_spike_flags);
        s1off_alignment_flags = ...
            valid_time >= s1off_alignment_onset - t1(s1_spike_flags) & ...
            valid_time < s1off_alignment_onset + max(xlim(sps(2)));
        s1off_chunk_flags = ...
            valid_time >= s1off_alignment_onset + min(xlim(sps(2))) & ...
            valid_time < s1off_alignment_onset + max(xlim(sps(2)));
        s1off_spkrates = s1_spike_rates;
        s1off_spkrates(~s1off_alignment_flags') = nan;
        s1off_spkrates = reshape(...
            s1off_spkrates(s1off_chunk_flags'),[s1off_n_tbins,s1_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(2))) + psthbin : psthbin : max(xlim(sps(2)));
        time_flags = time2plot <= max(xlim(sps(2)));
        
        % compute mean spike density function
        s1off_mu = nanmean(s1off_spkrates(:,time_flags),1);
        s1off_std = nanstd(s1off_spkrates(:,time_flags),0,1);
        s1off_trials_throughtime = sum(~isnan(s1off_spkrates),1);
        s1off_sem = s1off_std ./ sqrt(s1off_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(s1off_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        s1off_xpatch = [flagged_time,fliplr(flagged_time)];
        s1off_ypatch = [s1off_mu(~nan_flags)-s1off_sem(~nan_flags),...
            fliplr(s1off_mu(~nan_flags)+s1off_sem(~nan_flags))];
        patch_alpha = [s1off_trials_throughtime(~nan_flags),...
            fliplr(s1off_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(s1off_trials_throughtime)) / ...
            range(s1off_trials_throughtime);
        patch_alpha = patch_alpha * .25 + .05;
        patch(sps(2),s1off_xpatch,s1off_ypatch,i1_clrs(ii,:),...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'facealpha','interp',...
            'edgecolor','none');
        
        % patch average activity
        s1off_xpatch = flagged_time;
        s1off_ypatch = s1off_mu(~nan_flags);
        s1off_ypatch(end) = nan;
        patch_alpha = s1off_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(s1off_trials_throughtime)) / ...
            range(s1off_trials_throughtime);
        patch_alpha = patch_alpha * .95 + .05;
        patch(sps(2),s1off_xpatch,s1off_ypatch,0,...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'edgealpha','interp',...
            'edgecolor',i1_clrs(ii,:),...
            'linewidth',1.5);
        
        % plot stimulus offset
        plot(sps(2),time2plot(onset_flags),s1off_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',7.5,...
            'markerfacecolor',i1_clrs(ii,:),...
            'markeredgecolor','w');
        
        % plot I1 raster
        s1off_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(s1_spike_flags) + ...
            t1(s1_spike_flags));
        s1off_trial_idcs = (1 : s1_n_trials)' + s1_n_trial_counter;
        s1off_trial_mat = repmat(s1off_trial_idcs,1,n_paddedtimebins);
        s1off_spike_trials = s1off_trial_mat(s1_spike_counts >= 1);
        s1off_spike_times = s1off_time_mat(s1_spike_counts >= 1);
        trial_sorter = [t1(s1_spike_flags),prev_choice(s1_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + s1_n_trial_counter;
        s1off_sorted_trials = resorted_idcs(s1off_spike_trials - s1_n_trial_counter);
        plot(sps(2+n_cols),s1off_spike_times,s1off_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(2))) + [0,1,1,0] .* range(xlim(sps(2))) * .1/3;
        ypatch = [.5,.5,s1_n_trials+.5,s1_n_trials+.5] + s1_n_trial_counter;
        patch(sps(2+n_cols),xpatch,ypatch,i1_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % update trial counters
        s1_n_trial_counter = s1_n_trial_counter + s1_n_trials;
        s1_boundaries(ii) = s1_n_trials;
    end
    
    % iterate through intensities
    for ii = 1 : n_i
        i2_flags = i2 == i_set(ii);
        s2_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i2_flags;
        if sum(s2_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        s2_spike_counts = data.FR(s2_spike_flags,:);
        s2_spike_rates = conv2(...
            1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
        s2_n_trials = size(s2_spike_counts,1);
        s2on_n_tbins = range(xlim(sps(3)))  * psthbin;
        s2off_n_tbins = range(xlim(sps(4)))  * psthbin;
        
        % S2-onset-aligned, I2-split spike rates
        s2on_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            isi;
        s2on_alignment_flags = ...
            valid_time >= s2on_alignment_onset + min(xlim(sps(3))) & ...
            valid_time < s2on_alignment_onset + t2(s2_spike_flags);
        s2on_chunk_flags = ...
            valid_time >= s2on_alignment_onset + min(xlim(sps(3))) & ...
            valid_time < s2on_alignment_onset + max(xlim(sps(3)));
        s2on_spkrates = s2_spike_rates;
        s2on_spkrates(~s2on_alignment_flags') = nan;
        s2on_spkrates = reshape(...
            s2on_spkrates(s2on_chunk_flags'),[s2on_n_tbins,s2_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(3))) + psthbin : psthbin : max(xlim(sps(3)));
        time_flags = time2plot <= max(xlim(sps(3)));
        
        % compute mean spike density function
        s2on_mu = nanmean(s2on_spkrates(:,time_flags),1);
        s2on_std = nanstd(s2on_spkrates(:,time_flags),0,1);
        s2on_trials_throughtime = sum(~isnan(s2on_spkrates),1);
        s2on_sem = s2on_std ./ sqrt(s2on_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(s2on_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        s2on_xpatch = [flagged_time,fliplr(flagged_time)];
        s2on_ypatch = [s2on_mu(~nan_flags)-s2on_sem(~nan_flags),...
            fliplr(s2on_mu(~nan_flags)+s2on_sem(~nan_flags))];
        patch_alpha = [s2on_trials_throughtime(~nan_flags),...
            fliplr(s2on_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(s2on_trials_throughtime)) / ...
            range(s2on_trials_throughtime);
        patch_alpha = patch_alpha * .25 + .05;
        patch(sps(3),s2on_xpatch,s2on_ypatch,i2_clrs(ii,:),...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'facealpha','interp',...
            'edgecolor','none');
        
        % patch average activity
        s2on_xpatch = flagged_time;
        s2on_ypatch = s2on_mu(~nan_flags);
        s2on_ypatch(end) = nan;
        patch_alpha = s2on_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(s2on_trials_throughtime)) / ...
            range(s2on_trials_throughtime);
        patch_alpha = patch_alpha * .95 + .05;
        patch(sps(3),s2on_xpatch,s2on_ypatch,0,...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'edgealpha','interp',...
            'edgecolor',i2_clrs(ii,:),...
            'linewidth',1.5);
        
        % plot stimulus onset
        plot(sps(3),time2plot(onset_flags),s2on_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',7.5,...
            'markerfacecolor','w',...
            'markeredgecolor',i2_clrs(ii,:));
        
        % plot I2 raster
        s2on_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            isi);
        s2on_trial_idcs = (1 : s2_n_trials)' + s2_n_trial_counter;
        s2on_trial_mat = repmat(s2on_trial_idcs,1,n_paddedtimebins);
        s2on_spike_trials = s2on_trial_mat(s2_spike_counts >= 1);
        s2on_spike_times = s2on_time_mat(s2_spike_counts >= 1);
        trial_sorter = [t2(s2_spike_flags),choice(s2_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + s2_n_trial_counter;
        s2on_sorted_trials = resorted_idcs(s2on_spike_trials - s2_n_trial_counter);
        plot(sps(3+n_cols),s2on_spike_times,s2on_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(3))) + [0,1,1,0] .* range(xlim(sps(3))) * .1/3;
        ypatch = [.5,.5,s2_n_trials+.5,s2_n_trials+.5] + s2_n_trial_counter;
        patch(sps(3+n_cols),xpatch,ypatch,i2_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % S2-onset-aligned, I2-split spike rates
        s2off_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            isi + ...
            t2(s2_spike_flags);
        s2off_alignment_flags = ...
            valid_time >= s2off_alignment_onset - t2(s2_spike_flags) & ...
            valid_time < s2off_alignment_onset + max(xlim(sps(4)));
        s2off_chunk_flags = ...
            valid_time >= s2off_alignment_onset + min(xlim(sps(4))) & ...
            valid_time < s2off_alignment_onset + max(xlim(sps(4)));
        s2off_spkrates = s2_spike_rates;
        s2off_spkrates(~s2off_alignment_flags') = nan;
        s2off_spkrates = reshape(...
            s2off_spkrates(s2off_chunk_flags'),[s2off_n_tbins,s2_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(4))) + psthbin : psthbin : max(xlim(sps(4)));
        time_flags = time2plot <= max(xlim(sps(4)));
        
        % compute mean spike density function
        s2off_mu = nanmean(s2off_spkrates(:,time_flags),1);
        s2off_std = nanstd(s2off_spkrates(:,time_flags),0,1);
        s2off_trials_throughtime = sum(~isnan(s2off_spkrates),1);
        s2off_sem = s2off_std ./ sqrt(s2off_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(s2off_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        s2off_xpatch = [flagged_time,fliplr(flagged_time)];
        s2off_ypatch = [s2off_mu(~nan_flags)-s2off_sem(~nan_flags),...
            fliplr(s2off_mu(~nan_flags)+s2off_sem(~nan_flags))];
        patch_alpha = [s2off_trials_throughtime(~nan_flags),...
            fliplr(s2off_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(s2off_trials_throughtime)) / ...
            range(s2off_trials_throughtime);
        patch_alpha = patch_alpha * .25 + .05;
        patch(sps(4),s2off_xpatch,s2off_ypatch,i2_clrs(ii,:),...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'facealpha','interp',...
            'edgecolor','none');
        
        % patch average activity
        s2off_xpatch = flagged_time;
        s2off_ypatch = s2off_mu(~nan_flags);
        s2off_ypatch(end) = nan;
        patch_alpha = s2off_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(s2off_trials_throughtime)) / ...
            range(s2off_trials_throughtime);
        patch_alpha = patch_alpha * .95 + .05;
        patch(sps(4),s2off_xpatch,s2off_ypatch,0,...
            'facevertexalphadata',patch_alpha',...
            'alphadatamapping','none',...
            'edgealpha','interp',...
            'edgecolor',i2_clrs(ii,:),...
            'linewidth',1.5);
        
        % plot stimulus offset
        plot(sps(4),time2plot(onset_flags),s2off_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',7.5,...
            'markerfacecolor',i2_clrs(ii,:),...
            'markeredgecolor','w');
        
        % plot I2 raster
        s2off_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            isi + ...
            t2(s2_spike_flags));
        s2off_trial_idcs = (1 : s2_n_trials)' + s2_n_trial_counter;
        s2off_trial_mat = repmat(s2off_trial_idcs,1,n_paddedtimebins);
        s2off_spike_trials = s2off_trial_mat(s2_spike_counts >= 1);
        s2off_spike_times = s2off_time_mat(s2_spike_counts >= 1);
        trial_sorter = [t2(s2_spike_flags),choice(s2_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + s2_n_trial_counter;
        s2off_sorted_trials = resorted_idcs(s2off_spike_trials - s2_n_trial_counter);
        plot(sps(4+n_cols),s2off_spike_times,s2off_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(4))) + [0,1,1,0] .* range(xlim(sps(4))) * .1/3;
        ypatch = [.5,.5,s2_n_trials+.5,s2_n_trials+.5] + s2_n_trial_counter;
        patch(sps(4+n_cols),xpatch,ypatch,i2_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');
        
        % update trial counters
        s2_n_trial_counter = s2_n_trial_counter + s2_n_trials;
        s2_boundaries(ii) = s2_n_trials;
    end
    
    % text annotations
    %     mean_formula = strjoin(repmat({'%.1f'},1,n_i),' | ');
    %     text(sps(2),.05,1.3,sprintf(['mean FR: ',mean_formula],...
    %         nanmean(mean_frs(neurons2plot(nn),:,:),2)),...
    %         'horizontalalignment','left',...
    %         'verticalalignment','top',...
    %         'units','normalized');
    %     trialcount_formula = strjoin(repmat({'%i'},1,n_i),' | ');
    %     text(sps(2),.05,1.2,sprintf(['min. trial counts: ',trialcount_formula],...
    %         min(trial_type_counts(neurons2plot(nn),end-1,:),[],2)),...
    %         'horizontalalignment','left',...
    %         'verticalalignment','top',...
    %         'units','normalized');
    %     text(sps(2),.05,1.1,sprintf('stability coeff.: %.2f',...
    %         stability_coeffs(neurons2plot(nn))),...
    %         'horizontalalignment','left',...
    %         'verticalalignment','top',...
    %         'units','normalized');
    %     text(sps(2),1,1.3,sprintf('trial count flag: %i',...
    %         trial_count_flags(neurons2plot(nn))),...
    %         'horizontalalignment','right',...
    %         'verticalalignment','top',...
    %         'units','normalized');
    %     text(sps(2),1,1.2,sprintf('FR flag: %i',...
    %         mean_fr_flags(neurons2plot(nn))),...
    %         'horizontalalignment','right',...
    %         'verticalalignment','top',...
    %         'units','normalized');
    %     text(sps(2),1,1.1,sprintf('stability flag: %i',...
    %         stability_flags(neurons2plot(nn))),...
    %         'horizontalalignment','right',...
    %         'verticalalignment','top',...
    %         'units','normalized');
    
    if all(isnan(s1_boundaries)) || ...
            all(isnan(s2_boundaries))
        continue;
    end
    
    % axes updates
    set(sps([1,2]+n_cols),...
        'ylim',[1,nansum(s1_boundaries)]+[-1,1]*.5,...
        'ytick',unique([1;cumsum(s1_boundaries,'omitnan')]));
    set(sps([3,4]+n_cols),...
        'ylim',[1,nansum(s2_boundaries)]+[-1,1]*.5,...
        'ytick',unique([1;cumsum(s2_boundaries,'omitnan')]));
    
    % axes linkage
    %     linkaxes(sps(1:n_cols),'y');
    neurons2plot = [...
        38,68,72,215,234,381,391,393,402,436,...
        441,448,459,462,470,473,506,526];
    
    % update y-axis limits
    yylim_s1 = [];
    yylim_s2 = [];
    if ismember(neurons2plot(nn),[38])
        yylim_s1 = [5,45];
        yylim_s2 = [10,50];
    elseif ismember(neurons2plot(nn),[68,72])
        yylim_s1 = [0,15];
        yylim_s2 = [0,15];
    elseif ismember(neurons2plot(nn),[215])
        yylim_s1 = [0,65];
        yylim_s2 = [0,65];
    elseif ismember(neurons2plot(nn),[234])
        yylim_s1 = [10,50];
        yylim_s2 = [10,50];
    elseif ismember(neurons2plot(nn),[402])
        yylim_s1 = [0,45];
        yylim_s2 = [0,45];
    elseif ismember(neurons2plot(nn),[448])
        yylim_s1 = [10,45];
        yylim_s2 = [0,35];
    elseif ismember(neurons2plot(nn),[459])
        yylim_s1 = [0,50];
        yylim_s2 = [0,50];
    elseif ismember(neurons2plot(nn),[381])
        yylim_s1 = [5,40];
        yylim_s2 = [5,40];
    elseif ismember(neurons2plot(nn),[391])
        yylim_s1 = [0,20];
        yylim_s2 = [0,20];
    elseif ismember(neurons2plot(nn),[393])
        yylim_s1 = [15,45];
        yylim_s2 = [0,30];
    elseif ismember(neurons2plot(nn),[436])
        yylim_s1 = [0,30];
        yylim_s2 = [0,30];
    elseif ismember(neurons2plot(nn),[441])
        yylim_s1 = [0,25];
        yylim_s2 = [0,25];
    elseif ismember(neurons2plot(nn),[462])
        yylim_s1 = [25,60];
        yylim_s2 = [10,45];
    elseif ismember(neurons2plot(nn),[470])
        yylim_s1 = [0,15];
        yylim_s2 = [0,15];
    elseif ismember(neurons2plot(nn),[473])
        yylim_s1 = [25,50];
        yylim_s2 = [0,25];
    elseif ismember(neurons2plot(nn),[506])
        yylim_s1 = [0,30];
        yylim_s2 = [0,30];
    elseif ismember(neurons2plot(nn),[526])
        yylim_s1 = [0,50];
        yylim_s2 = [0,50];
    end
    if ~isempty(yylim_s1)
        set(sps(1:n_cols/2),...
            'ylim',yylim_s1,...
            'ytick',yylim_s1,...
            'clipping','off');
    end
    if ~isempty(yylim_s2)
        set(sps(n_cols/2+1:n_cols),...
            'ylim',yylim_s2,...
            'ytick',yylim_s2,...
            'clipping','off');
    end
    drawnow;
    
    % update y-labels
    ylbl_5.Position(1) = ylbl_1.Position(1);
    ylbl_6.Position(1) = ylbl_2.Position(1);
    ylbl_7.Position(1) = ylbl_3.Position(1);
    ylbl_8.Position(1) = ylbl_4.Position(1);
    
    % shift psth axes down a bit
    for ii = 1 : n_cols
        set(sps(ii),...
            'position',get(sps(ii),'position')-[0,.075,0,0]);
    end
    
    % plot waveform
    %     plot(axwav,data.Shape(neurons2plot(nn),:),...
    %         'color','k',...
    %         'linewidth',1.5);
    
    % plot whether or not this passed selection
    %     if ismember(neurons2plot(nn),flagged_neurons)
    %         clr = [0,1,0];
    %     else
    %         clr = [1,0,0];
    %     end
    %     plot(axwav,max(xlim(axwav)),min(ylim(axwav)),...
    %         'marker','.',...
    %         'markersize',25,...
    %         'color',clr);
    
    % save figure
    if want2save
        try
            % save settings
            png_file = fullfile(raster_path,[get(fig,'name'),'.png']);
            print(fig,png_file,'-dpng','-r600','-painters');
            %             svg_file = fullfile(panel_path,[fig.Name,'.svg']);
            %             print(fig,svg_file,'-dsvg','-painters');
            close(fig);
        catch
            close(fig);
        end
    else
        pause(1);
        close(fig);
    end
end
