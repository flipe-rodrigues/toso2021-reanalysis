%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% neuron selection

% selected for being good examples of i2-modulation
if strcmpi(task_str,'duration')
    neurons2plot = fliplr([...
        393,215,72,526]);
elseif strcmpi(task_str,'intensity')
    neurons2plot = [...
        19,22,30,61,66,70,100,111,112,115,...
        166,238,243,260,344,408,410];
end
n_neurons2plot = numel(neurons2plot);

%% ROI settings

% preallocation
rois = struct();
n_bins = struct();
time = struct();
psths = struct();

% roi definition
rois.s2on = [0,t_set(end)];
rois.s2off = sort(-rois.s2on);

% iterate through task epochs
task_epochs = fieldnames(rois);
n_epochs = numel(task_epochs);
for ii = 1 : n_epochs
    epoch = task_epochs{ii};
    n_bins.(epoch) = range(rois.(epoch)) / psthbin;
    time.(epoch) = linspace(rois.(epoch)(1),rois.(epoch)(2),n_bins.(epoch));
    psths.(epoch) = nan(n_bins.(epoch),n_neurons_total,n_i);
end

%% construct T1-aligned, Ii-split psths

% half-width smoothing kernel
hw_kernel = gammakernel('peakx',kernel_peak_time/2,'binwidth',psthbin);

% iterate through neurons
for nn = 1 : n_neurons2plot
    progressreport(nn,n_neurons2plot,'fetching middle-I2 SDFs');
    neuron_flags = data.NeuronNumb == neurons2plot(nn);
    
    % trial selection
    i2_flags = i2 == i_set(i2_mode_idx);
    trial_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i2_flags;
    if sum(trial_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    spike_counts = data.FR(trial_flags,:);
    spike_rates = conv2(...
        1,hw_kernel.pdf,spike_counts / psthbin * 1e3,'same');
    spike_rates = conv2(...
        1,fliplr(hw_kernel.pdf),spike_rates,'same');
    n_trials = size(spike_counts,1);
    
    % S2-onset-aligned spike rates
    s2on_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(trial_flags) + ...
        t1(trial_flags) + ...
        isi;
    s2on_alignment_flags = ...
        padded_time >= s2on_alignment + rois.s2on(1) & ...
        padded_time < s2on_alignment + t2(trial_flags);
    s2on_chunk_flags = ...
        padded_time >= s2on_alignment + rois.s2on(1) & ...
        padded_time < s2on_alignment + rois.s2on(2);
    s2on_spkrates = spike_rates';
    s2on_spkrates(~s2on_alignment_flags') = nan;
    s2on_spkrates = reshape(...
        s2on_spkrates(s2on_chunk_flags'),[n_bins.s2on,n_trials])';
    
    % S2-offset-aligned spike rates
    s2off_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(trial_flags) + ...
        t1(trial_flags) + ...
        isi + ...
        t2(trial_flags);
    s2off_alignment_flags = ...
        padded_time >= s2off_alignment - t2(trial_flags) & ...
        padded_time < s2off_alignment + rois.s2off(2);
    s2off_chunk_flags = ...
        padded_time >= s2off_alignment + rois.s2off(1) & ...
        padded_time < s2off_alignment + rois.s2off(2);
    s2off_spkrates = spike_rates';
    s2off_spkrates(~s2off_alignment_flags') = nan;
    s2off_spkrates = reshape(...
        s2off_spkrates(s2off_chunk_flags'),[n_bins.s2off,n_trials])';
    
    % compute mean spike density function
    psths.s2on(:,nn,i2_mode_idx) = nanmean(s2on_spkrates,1);
    psths.s2off(:,nn,i2_mode_idx) = nanmean(s2off_spkrates,1);
end

%% generate I2-modulated spike rates

% modulation settings
modulation = log(i_set) ./ log(i_set(i2_mode_idx));
scaling = modulation .^ 0;
gain = modulation .^ 1 * 1;
offset = modulation .^ 1 * 0;

% iterate through neurons
for nn = 1 : n_neurons2plot
    progressreport(nn,n_neurons2plot,'simulating I2 modulation');

    % iterate through intensities
    for ii = 1 : n_i
        if ii == i2_mode_idx
            continue;
        end

        % apply I2 modulation at S2 presentation
        psths.s2on(:,nn,ii) = offset(ii) + gain(ii) * ...
            interp1(time.s2on,psths.s2on(:,nn,i2_mode_idx),time.s2on*scaling(ii),...
            'linear','extrap');
        psths.s2off(:,nn,ii) = offset(ii) + gain(ii) * ...
            interp1(time.s2off,psths.s2off(:,nn,i2_mode_idx),time.s2off*scaling(ii));
    end
end

%% construct Si-aligned, Ti- & Ii-split psths
ti_padd = [-0,0];
sdf_gain = .75;
sdf_spacing = 1.05;

% figure initialization
fig = figure(figopt,...
    'position',[489,343,460,420],...
    'name',sprintf('fake_modulation_%s',contrast_str),...
    'color',[1,1,1]*245/255);

% axes initialization
sps = gobjects(2,1);
sps(1) = subplot(1,2,1);
sps(2) = subplot(1,2,2);
set(sps,...
    axesopt.default,...
    'plotboxaspectratiomode','auto',...
    'clipping','off',...
    'ylim',[0,n_neurons2plot+1]+[-1,1]*.05*(n_neurons2plot+1),...
    'ycolor','none');
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
    'xticklabel',xxticklabel);
xlabel(sps(1),'Time since S2 onset (ms)');
xlabel(sps(2),'Time since S2 offset (ms)');
ylabel(sps(1),'Firing rate (Hz)');
ylabel(sps(2),'Firing rate (Hz)');

% iterate through neurons
for nn = 1 : n_neurons2plot
    progressreport(nn,n_neurons2plot,'plotting simulated SDFs');
    neuron_flags = data.NeuronNumb == neurons2plot(nn);

    % iterate through intensities
    for ii = n_i : -1 : 1
        i2_flags = i2 == i_set(ii);
        
        % trial selection
        trial_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i2_flags;
        n_trials = sum(trial_flags);
        if n_trials == 0
            continue;
        end
        
        % compute surviving trial counts (through time)
        time2plot = min(xlim(sps(1))) + psthbin : psthbin : max(xlim(sps(1)));
        time_mat = repmat(time2plot,n_total_trials,1);
        s2on_surviving_trial_counts = ...
            sum(time_mat(trial_flags,:) <= t2(trial_flags));

        % patch average activity
        s2on_mu = psths.s2on(:,nn,ii)';
        nan_flags = isnan(s2on_mu);
        if ii == n_i
            s2_min = min(s2on_mu);
            s2_range = range(s2on_mu);
        end
        s2on_mu = (s2on_mu - s2_min) / s2_range - .5;
        s2on_mu = s2on_mu * sdf_gain + .5 + (nn - 1) * sdf_spacing;
        s2on_xpatch = time2plot(~nan_flags);
        s2on_ypatch = s2on_mu(~nan_flags);
        s2on_ypatch(end) = nan;
        patch_alpha = s2on_surviving_trial_counts(~nan_flags);
        patch_alpha = patch_alpha ./ max(patch_alpha) .* ...
            range(alphabounds_mu) + alphabounds_mu(1);
        alpha_levels = unique(patch_alpha,'stable');
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(1),...
                [s2on_xpatch(alpha_flags),nan],...
                [s2on_ypatch(alpha_flags),nan],0,...
                'edgealpha',alpha_levels(aa),...
                'edgecolor',i2_clrs(ii,:),...
                'facecolor','none',...
                'linewidth',1.5);
        end
        
        % compute surviving trial counts (through time)
        time2plot = min(xlim(sps(2))) + psthbin : psthbin : max(xlim(sps(2)));
        time_mat = repmat(time2plot,n_total_trials,1);
        s2off_surviving_trial_counts = ...
            sum(time_mat(trial_flags,:) > -t2(trial_flags));
        
        % patch average activity
        s2off_mu = psths.s2off(:,nn,ii)';
        nan_flags = isnan(s2off_mu);
        s2off_mu = (s2off_mu - s2_min) / s2_range - .5;
        s2off_mu = s2off_mu * sdf_gain + .5 + (nn - 1) * sdf_spacing;
        s2off_xpatch = time2plot(~nan_flags);
        s2off_ypatch = s2off_mu(~nan_flags);
        s2off_ypatch(end) = nan;
        patch_alpha = s2off_surviving_trial_counts(~nan_flags);
        patch_alpha = patch_alpha ./ max(patch_alpha) .* ...
            range(alphabounds_mu) + alphabounds_mu(1);
        alpha_levels = unique(patch_alpha,'stable');
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(2),...
                [s2off_xpatch(alpha_flags),nan],...
                [s2off_ypatch(alpha_flags),nan],0,...
                'edgealpha',alpha_levels(aa),...
                'edgecolor',i2_clrs(ii,:),...
                'facecolor','none',...
                'linewidth',1.5);
        end
    end
end

% iterate through durations
for tt = 1 : n_t
    multiplier = 1 - .015 * (tt - 1);
    
    % plot spike integration window
    plot(sps(1),[0,t_set(tt)],[1,1]*max(ylim(sps(1)))*multiplier,'-k',...
        'linewidth',3);
    plot(sps(2),[-t_set(tt),0],[1,1]*max(ylim(sps(2)))*multiplier,'-k',...
        'linewidth',3);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end