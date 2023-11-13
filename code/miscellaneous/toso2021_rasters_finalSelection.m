%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end
close all;

%% neuron selection

% selected for being good examples of i2-modulation
if strcmpi(task_str,'duration')
    neurons2plot = fliplr([...
        393,473,215,72,526]);
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

% figure initialization
fig = figure(figopt,...
    ...'windowstate','maximized',...
    'position',[725,905,2100,1317],...
    'name','rasters_egcells');

% axes initialization;
n_rows = 4;
n_cols = n_neurons2plot;
n_sps = n_rows * n_cols;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(n_rows,n_cols,ii);
end
xxtick = unique([ti_padd(1);-pre_s1_delay;0;t_set;t_set(end)+ti_padd(2)]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < 1e3) = {''};
set(sps,...
    axesopt.default,...
    'layer','top',...
    'plotboxaspectratiomode','auto',...
    'xlim',[0,max(t_set)]+ti_padd,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight');
set(sps((1:n_cols)+[0;2]*n_cols),...
    'plotboxaspectratio',[2,1,1],...
    'xcolor','none');
for ii = 1 : n_cols
    xlabel(sps(ii+n_cols*1),'Time since S_1 onset (s)');
    xlabel(sps(ii+n_cols*3),'Time since S_2 onset (s)');
    ylabel(sps(ii+n_cols*0),'Firing rate (Hz)');
    ylabel(sps(ii+n_cols*2),'Firing rate (Hz)');
    ylabel(sps(ii+n_cols*1),'Trial #');
    ylabel(sps(ii+n_cols*3),'Trial #');
end

% clamping
i1_clamp_flags = i1 == i_set(i1_mode_idx);
i2_clamp_flags = i2 == i_set(i2_mode_idx);

% iterate through neurons
for nn = 1 : n_neurons2plot
    progressreport(nn,n_neurons2plot,'parsing neural data');
    neuron_flags = data.NeuronNumb == neurons2plot(nn);

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
        s1_n_tbins = range(xlim(sps(nn+n_cols*0))) * psthbin;
        
        % S1-onset-aligned, I1-split spike rates
        s1_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(s1_spike_flags);
        s1_alignment_flags = ...
            valid_time >= s1_alignment_onset + min(xlim(sps(nn+n_cols*0))) & ...
            valid_time < s1_alignment_onset + t1(s1_spike_flags);
        s1_chunk_flags = ...
            valid_time >= s1_alignment_onset + min(xlim(sps(nn+n_cols*0))) & ...
            valid_time < s1_alignment_onset + max(xlim(sps(nn+n_cols*0)));
        s1_spkrates = s1_spike_rates;
        s1_spkrates(~s1_alignment_flags') = nan;
        s1_spkrates = reshape(...
            s1_spkrates(s1_chunk_flags'),[s1_n_tbins,s1_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(nn+n_cols*0))) + psthbin : psthbin : max(xlim(sps(nn+n_cols*0)));
        time_flags = time2plot <= max(xlim(sps(nn+n_cols*0)));
        
        % compute mean spike density function
        s1_mu = nanmean(s1_spkrates(:,time_flags),1);
        s1_std = nanstd(s1_spkrates(:,time_flags),0,1);
        s1_trials_throughtime = sum(~isnan(s1_spkrates),1);
        s1_sem = s1_std ./ sqrt(s1_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(s1_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        s1_xpatch = [flagged_time,fliplr(flagged_time)];
        s1_ypatch = [s1_mu(~nan_flags)-s1_sem(~nan_flags),...
            fliplr(s1_mu(~nan_flags)+s1_sem(~nan_flags))];
        patch_alpha = [s1_trials_throughtime(~nan_flags),...
            fliplr(s1_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(s1_trials_throughtime)) / ...
            range(s1_trials_throughtime);
        patch_alpha = patch_alpha * .25 + .05;
        alpha_levels = unique(patch_alpha);
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(nn+n_cols*0),...
                s1_xpatch(alpha_flags),...
                s1_ypatch(alpha_flags),0,...
                'facealpha',alpha_levels(aa),...
                'edgecolor','none',...
                'facecolor',i1_clrs(ii,:));
        end
        
        % patch average activity
        s1_xpatch = flagged_time;
        s1_ypatch = s1_mu(~nan_flags);
        s1_ypatch(end) = nan;
        patch_alpha = s1_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(s1_trials_throughtime)) / ...
            range(s1_trials_throughtime);
        patch_alpha = patch_alpha * .95 + .05;
        alpha_levels = unique(patch_alpha);
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(nn+n_cols*0),...
                [s1_xpatch(alpha_flags),nan],...
                [s1_ypatch(alpha_flags),nan],0,...
                'edgealpha',alpha_levels(aa),...
                'edgecolor',i1_clrs(ii,:),...
                'facecolor','none',...
                'linewidth',1.5);
        end
        
        % plot stimulus onset
        plot(sps(nn+n_cols*0),time2plot(onset_flags),s1_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',7.5,...
            'markerfacecolor','w',...
            'markeredgecolor',i1_clrs(ii,:));
        
        % plot I1 raster
        s1_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(s1_spike_flags));
        s1_trial_idcs = (1 : s1_n_trials)' + s1_n_trial_counter;
        s1_trial_mat = repmat(s1_trial_idcs,1,n_paddedtimebins);
        s1_spike_trials = s1_trial_mat(s1_spike_counts >= 1);
        s1_spike_times = s1_time_mat(s1_spike_counts >= 1);
        trial_sorter = [t1(s1_spike_flags),prev_choice(s1_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + s1_n_trial_counter;
        s1_sorted_trials = resorted_idcs(s1_spike_trials - s1_n_trial_counter);
        plot(sps(nn+n_cols*1),s1_spike_times,s1_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(nn+n_cols*0))) + [0,1,1,0] .* range(xlim(sps(nn+n_cols*0))) * .1/3;
        ypatch = [.5,.5,s1_n_trials+.5,s1_n_trials+.5] + s1_n_trial_counter;
        patch(sps(nn+n_cols*1),xpatch,ypatch,i1_clrs(ii,:),...
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
        s2_n_tbins = range(xlim(sps(nn+n_cols*2)))  * psthbin;
        
        % S2-onset-aligned, I2-split spike rates
        s2_alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            isi;
        s2_alignment_flags = ...
            valid_time >= s2_alignment_onset + min(xlim(sps(nn+n_cols*2))) & ...
            valid_time < s2_alignment_onset + t2(s2_spike_flags);
        s2_chunk_flags = ...
            valid_time >= s2_alignment_onset + min(xlim(sps(nn+n_cols*2))) & ...
            valid_time < s2_alignment_onset + max(xlim(sps(nn+n_cols*2)));
        s2_spkrates = s2_spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(...
            s2_spkrates(s2_chunk_flags'),[s2_n_tbins,s2_n_trials])';
        
        % time selection
        time2plot = min(xlim(sps(nn+n_cols*2))) + psthbin : psthbin : max(xlim(sps(nn+n_cols*2)));
        time_flags = time2plot <= max(xlim(sps(nn+n_cols*2)));
        
        % compute mean spike density function
        s2_mu = nanmean(s2_spkrates(:,time_flags),1);
        s2_std = nanstd(s2_spkrates(:,time_flags),0,1);
        s2_trials_throughtime = sum(~isnan(s2_spkrates),1);
        s2_sem = s2_std ./ sqrt(s2_trials_throughtime);
        
        % flag current stimulus period
        nan_flags = isnan(s2_mu);
        onset_flags = time2plot <= 0 & [time2plot(2:end),nan] > 0;
        offset_flags = diff([nan_flags,true]) == 1;
        flagged_time = time2plot(~nan_flags);
        
        % patch s.e.m.
        s2_xpatch = [flagged_time,fliplr(flagged_time)];
        s2_ypatch = [s2_mu(~nan_flags)-s2_sem(~nan_flags),...
            fliplr(s2_mu(~nan_flags)+s2_sem(~nan_flags))];
        patch_alpha = [s2_trials_throughtime(~nan_flags),...
            fliplr(s2_trials_throughtime(~nan_flags))];
        patch_alpha = (patch_alpha - min(s2_trials_throughtime)) / ...
            range(s2_trials_throughtime);
        patch_alpha = patch_alpha * .25 + .05;
        alpha_levels = unique(patch_alpha);
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(nn+n_cols*2),...
                s2_xpatch(alpha_flags),...
                s2_ypatch(alpha_flags),0,...
                'facealpha',alpha_levels(aa),...
                'edgecolor','none',...
                'facecolor',i2_clrs(ii,:));
        end
        
        % patch average activity
        s2_xpatch = flagged_time;
        s2_ypatch = s2_mu(~nan_flags);
        s2_ypatch(end) = nan;
        patch_alpha = s2_trials_throughtime(~nan_flags);
        patch_alpha = (patch_alpha - min(s2_trials_throughtime)) / ...
            range(s2_trials_throughtime);
        patch_alpha = patch_alpha * .95 + .05;
        alpha_levels = unique(patch_alpha);
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(nn+n_cols*2),...
                [s2_xpatch(alpha_flags),nan],...
                [s2_ypatch(alpha_flags),nan],0,...
                'edgealpha',alpha_levels(aa),...
                'edgecolor',i2_clrs(ii,:),...
                'facecolor','none',...
                'linewidth',1.5);
        end
        
        % plot stimulus onset
        plot(sps(nn+n_cols*2),time2plot(onset_flags),s2_mu(onset_flags),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',7.5,...
            'markerfacecolor','w',...
            'markeredgecolor',i2_clrs(ii,:));
        
        % plot I2 raster
        s2_time_mat = padded_time - (...
            pre_init_padding + ...
            pre_s1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            isi);
        s2_trial_idcs = (1 : s2_n_trials)' + s2_n_trial_counter;
        s2_trial_mat = repmat(s2_trial_idcs,1,n_paddedtimebins);
        s2_spike_trials = s2_trial_mat(s2_spike_counts >= 1);
        s2_spike_times = s2_time_mat(s2_spike_counts >= 1);
        trial_sorter = [t2(s2_spike_flags),choice(s2_spike_flags)];
        [~,sorted_idcs] = sortrows(trial_sorter,[1,-2]);
        [~,resorted_idcs] = sortrows(sorted_idcs);
        resorted_idcs = resorted_idcs + s2_n_trial_counter;
        s2_sorted_trials = resorted_idcs(s2_spike_trials - s2_n_trial_counter);
        plot(sps(nn+n_cols*3),s2_spike_times,s2_sorted_trials,...
            'color','k',...
            'marker',spike_marker,...
            'markersize',spike_markersize,...
            'linestyle','none');
        
        % plot raster bands
        xpatch = min(xlim(sps(nn+n_cols*2))) + [0,1,1,0] .* range(xlim(sps(nn+n_cols*2))) * .1/3;
        ypatch = [.5,.5,s2_n_trials+.5,s2_n_trials+.5] + s2_n_trial_counter;
        patch(sps(nn+n_cols*3),xpatch,ypatch,i2_clrs(ii,:),...
            'linewidth',1.5,...
            'facealpha',.75,...
            'edgecolor','none');

        % update trial counters
        s2_n_trial_counter = s2_n_trial_counter + s2_n_trials;
        s2_boundaries(ii) = s2_n_trials;
    end

    % axes updates
    set(sps(nn+n_cols*1),...
        'ylim',[1,nansum(s1_boundaries)]+[-1,1]*.5,...
        'ytick',unique([1;cumsum(s1_boundaries,'omitnan')]));
    set(sps(nn+n_cols*3),...
        'ylim',[1,nansum(s2_boundaries)]+[-1,1]*.5,...
        'ytick',unique([1;cumsum(s2_boundaries,'omitnan')]));

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
        set(sps(nn+n_cols*0),...
            'ylim',yylim_s1,...
            'ytick',yylim_s1,...
            'clipping','off');
    end
    if ~isempty(yylim_s2)
        set(sps(nn+n_cols*2),...
            'ylim',yylim_s2,...
            'ytick',yylim_s2,...
            'clipping','off');
    end
    drawnow;
end

% make sure psth & raster axes are the same size
for ii = 1 : n_cols
    for jj = [0,2]
        pos_psth = get(sps(ii+n_cols*(0+jj)),'position');
        pos_raster = get(sps(ii+n_cols*(1+jj)),'position');
        set(sps(ii+n_cols*(0+jj)),...
            'position',[pos_raster(1),pos_psth(2),pos_raster(3),pos_raster(4)]);
    end
end

% shift psth axes down a bit
for ii = 1 : n_cols
    set(sps(ii+n_cols*0),...
        'position',get(sps(ii+n_cols*0),'position')-[0,.065,0,0]);
    set(sps(ii+n_cols*2),...
        'position',get(sps(ii+n_cols*2),'position')-[0,.065,0,0]);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end