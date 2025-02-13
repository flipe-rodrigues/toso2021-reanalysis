%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% ROI settings

% preallocation
rois = struct();
n_bins = struct();
time = struct();
psths = struct();

% smoothing settings
gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);
gauss_padded_time = ...
    (1 : psthbin : n_paddedtimebins * psthbin) - psthbin;
gauss_validtime_flags = ...
    gauss_padded_time >= gauss_padded_time(1) - gauss_kernel.paddx(1) & ...
    gauss_padded_time <= gauss_padded_time(end) - gauss_kernel.paddx(end) + psthbin;
gauss_valid_time = gauss_padded_time(gauss_validtime_flags);

% roi definition
rois.pre_s1 = [gauss_kernel.paddx(1),unique(pre_s1_delay(valid_flags))];
rois.s1 = [0,t_set(end)];
rois.isi = [0,isi];
rois.s2on = [0,t_set(end)];
rois.s2off = sort(-rois.s2on);
rois.post_s2 = [0,post_s2_delay];
rois.go = [0,gauss_kernel.paddx(2)];

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

% iterate through neurons
for nn = 1 : n_neurons_total
    progressreport(nn,n_neurons_total,'parsing neural data');
    
    % trial selection
    neuron_flags = data.NeuronNumb == nn;
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
        1,gauss_kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
    n_trials = size(spike_counts,1);
    
    % pre S1 delay-aligned spike rates
    pre_s1_alignment = ...
        repmat(pre_init_padding,n_trials,1);
    pre_s1_alignment_flags = ...
        gauss_valid_time >= pre_s1_alignment + rois.pre_s1(1) & ...
        gauss_valid_time < pre_s1_alignment + pre_s1_delay(spike_flags);
    pre_s1_chunk_flags = ...
        gauss_valid_time >= pre_s1_alignment + rois.pre_s1(1) & ...
        gauss_valid_time < pre_s1_alignment + rois.pre_s1(2);
    pre_s1_spkrates = spike_rates;
    pre_s1_spkrates(~pre_s1_alignment_flags') = nan;
    pre_s1_spkrates = reshape(...
        pre_s1_spkrates(pre_s1_chunk_flags'),[n_bins.pre_s1,n_trials])';
    
    % S1-aligned spike rates
    s1_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags);
    s1_alignment_flags = ...
        gauss_valid_time >= s1_alignment + rois.s1(1) & ...
        gauss_valid_time < s1_alignment + t1(spike_flags);
    s1_chunk_flags = ...
        gauss_valid_time >= s1_alignment + rois.s1(1) & ...
        gauss_valid_time < s1_alignment + rois.s1(2);
    s1_spkrates = spike_rates;
    s1_spkrates(~s1_alignment_flags') = nan;
    s1_spkrates = reshape(...
        s1_spkrates(s1_chunk_flags'),[n_bins.s1,n_trials])';
    
    % inter S1-S2 delay-aligned spike rates
    isi_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags);
    isi_alignment_flags = ...
        gauss_valid_time >= isi_alignment + rois.isi(1) & ...
        gauss_valid_time < isi_alignment + isi;
    isi_chunk_flags = ...
        gauss_valid_time >= isi_alignment + rois.isi(1) & ...
        gauss_valid_time < isi_alignment + rois.isi(2);
    isi_spkrates = spike_rates;
    isi_spkrates(~isi_alignment_flags') = nan;
    isi_spkrates = reshape(...
        isi_spkrates(isi_chunk_flags'),[n_bins.isi,n_trials])';
    
    % S2-onset-aligned spike rates
    s2on_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi;
    s2on_alignment_flags = ...
        gauss_valid_time >= s2on_alignment + rois.s2on(1) & ...
        gauss_valid_time < s2on_alignment + t2(spike_flags);
    s2on_chunk_flags = ...
        gauss_valid_time >= s2on_alignment + rois.s2on(1) & ...
        gauss_valid_time < s2on_alignment + rois.s2on(2);
    s2on_spkrates = spike_rates;
    s2on_spkrates(~s2on_alignment_flags') = nan;
    s2on_spkrates = reshape(...
        s2on_spkrates(s2on_chunk_flags'),[n_bins.s2on,n_trials])';
    
    % S2-offset-aligned spike rates
    s2off_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi + ...
        t2(spike_flags);
    s2off_alignment_flags = ...
        gauss_valid_time >= s2off_alignment - t2(spike_flags) & ...
        gauss_valid_time < s2off_alignment + rois.s2off(2);
    s2off_chunk_flags = ...
        gauss_valid_time >= s2off_alignment + rois.s2off(1) & ...
        gauss_valid_time < s2off_alignment + rois.s2off(2);
    s2off_spkrates = spike_rates;
    s2off_spkrates(~s2off_alignment_flags') = nan;
    s2off_spkrates = reshape(...
        s2off_spkrates(s2off_chunk_flags'),[n_bins.s2off,n_trials])';
    
    % post S2 delay-aligned spike rates
    post_s2_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi + ...
        t2(spike_flags);
    post_s2_alignment_flags = ...
        gauss_valid_time >= post_s2_alignment + rois.post_s2(1) & ...
        gauss_valid_time < post_s2_alignment + post_s2_delay;
    post_s2_chunk_flags = ...
        gauss_valid_time >= post_s2_alignment + rois.post_s2(1) & ...
        gauss_valid_time < post_s2_alignment + rois.post_s2(2);
    post_s2_spkrates = spike_rates;
    post_s2_spkrates(~post_s2_alignment_flags') = nan;
    post_s2_spkrates = reshape(...
        post_s2_spkrates(post_s2_chunk_flags'),[n_bins.post_s2,n_trials])';
    
    % go cue-aligned spike rates
    go_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi + ...
        t2(spike_flags) + ...
        post_s2_delay;
    go_alignment_flags = ...
        gauss_valid_time >= go_alignment + rois.go(1) & ...
        gauss_valid_time < go_alignment + rois.go(2);
    go_chunk_flags = go_alignment_flags;
    go_spkrates = spike_rates;
    go_spkrates(~go_alignment_flags') = nan;
    go_spkrates = reshape(...
        go_spkrates(go_chunk_flags'),[n_bins.go,n_trials])';
    
    % compute mean spike density function
    psths.pre_s1(:,nn,i2_mode_idx) = nanmean(pre_s1_spkrates,1);
    psths.s1(:,nn,i2_mode_idx) = nanmean(s2on_spkrates,1); % !!! IMPORTANT DETAIL: S2, NOT S1
    psths.isi(:,nn,i2_mode_idx) = nanmean(isi_spkrates,1);
    psths.s2on(:,nn,i2_mode_idx) = nanmean(s2on_spkrates,1);
    psths.s2off(:,nn,i2_mode_idx) = nanmean(s2off_spkrates,1);
    psths.post_s2(:,nn,i2_mode_idx) = nanmean(post_s2_spkrates,1);
    psths.go(:,nn,i2_mode_idx) = nanmean(go_spkrates,1);
end

%% generate I2-modulated spike rates

% modulation settings
modulation = log(i_set) ./ log(i_set(i2_mode_idx));
scaling = modulation .^ 0;
gain = modulation .^ 1 * 1;
offset = modulation .^ 1 * 0;

% iterate through neurons
for nn = 1 : n_neurons_total
    progressreport(nn,n_neurons_total,'generating fake rates');
    if all(isnan(psths.s2on(:,nn,i2_mode_idx)))
        continue;
    end
    
    % iterate through intensities
    for ii = 1 : n_i
        if ii == i2_mode_idx
            continue;
        end
        
        % apply I1 modulation at S1 presentation
        psths.s1(:,nn,ii) = 1 * ...
            psths.s1(:,nn,i1_mode_idx);
        
        % apply I2 modulation at S2 presentation
        psths.s2on(:,nn,ii) = offset(ii) + gain(ii) * ...
            interp1(time.s2on,psths.s2on(:,nn,i2_mode_idx),time.s2on*scaling(ii),...
            'linear','extrap');
        psths.s2off(:,nn,ii) = offset(ii) + gain(ii) * ...
            interp1(time.s2off,psths.s2off(:,nn,i2_mode_idx),time.s2off*scaling(ii));

        % no intensity modulation in the remaining epochs
        psths.pre_s1(:,nn,ii) = psths.pre_s1(:,nn,i2_mode_idx);
        psths.isi(:,nn,ii) = psths.isi(:,nn,i2_mode_idx);
        psths.post_s2(:,nn,ii) = psths.post_s2(:,nn,i2_mode_idx);
        psths.go(:,nn,ii) = psths.go(:,nn,i2_mode_idx);
    end
end

%% generate I2-modulated spike trains

% preallocation
if ~isfield(data,'FakeFR')
    data.FakeFR = nan(size(data.FR));
end

% iterate through neurons
for nn = 1 : n_neurons_total
    progressreport(nn,n_neurons_total,'generating fake spikes');
    neuron_flags = data.NeuronNumb == nn;
    
    % trial selection
    spike_flags = ...
        valid_flags & ...
        neuron_flags;
    spike_trials = find(spike_flags);
    n_trials = numel(spike_trials);
    
    if ismember(nn,flagged_neurons)
        figure;
        hold on;
    end
    
    % iterate through trials
    for kk = 1 : n_trials
        trial_idx = spike_trials(kk);

        lambda_pre_s1 = psths.pre_s1(:,nn,i1(trial_idx)==i_set);
        lambda_s1 = psths.s1(time.s1<=t1(trial_idx),nn,i1(trial_idx)==i_set) + ...
            lambda_pre_s1(end) - psths.s1(1,nn,i1(trial_idx)==i_set);
        lambda_isi = psths.isi(:,nn,i1(trial_idx)==i_set) + ...
            lambda_s1(end) - psths.isi(1,nn,i1(trial_idx)==i_set);
        lambda_s2 = psths.s2on(time.s2on<=t2(trial_idx),nn,i2(trial_idx)==i_set) + ...
            lambda_isi(end) - psths.s2on(1,nn,i2(trial_idx)==i_set);
        lambda_post_s2 = psths.post_s2(:,nn,i2(trial_idx)==i_set) + ...
            lambda_s2(end) - psths.post_s2(1,nn,i2(trial_idx)==i_set);
        lambda_go = psths.go(:,nn,i2(trial_idx)==i_set) + ...
            lambda_post_s2(end) - psths.go(1,nn,i2(trial_idx)==i_set);
        
        lambda = [...
            lambda_pre_s1;...
            lambda_s1;...
            lambda_isi;...
            lambda_s2;...
            lambda_post_s2;...
            lambda_go];
        
        dur = numel(lambda) * psthbin;
        
        if ismember(nn,flagged_neurons)
            plot(lambda,...
                'color',i2_clrs(i2(trial_idx)==i_set,:));
            a=1
        end
        
        [n,ts] = poissonprocess(lambda,dur / 1e3);
        offset = pre_init_padding + rois.pre_s1(1);
        data.FakeFR(trial_idx,:) = ...
            histcounts(ts*1e3,[gauss_padded_time,n_paddedtimebins*psthbin]-offset);
    end
end