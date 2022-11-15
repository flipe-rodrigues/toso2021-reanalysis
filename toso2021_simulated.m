%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% ROI settings

% preallocation
rois = struct();
n_bins = struct();
time = struct();
psths = struct();

% smoothing settings
gauss_kernel = gausskernel('mu',0,'sig',50,'binwidth',psthbin);
    
% roi definition
rois.pre_s1 = [gauss_kernel.paddx(1),unique(pre_s1_delay(valid_flags))];
rois.s1 = [0,t_set(end)];
rois.inter_s1s2 = [0,isi];
rois.s2 = [0,t_set(end)];
rois.post_s2 = [0,post_s2_delay];
rois.go = [0,gauss_kernel.paddx(2)];

% roi alignment labels
alignments.pre_s1 = 'pre-T_{1}-delay';
alignments.s1 = 'T_{1}';
alignments.inter_s1s2 = 'inter T_{1}-T_{2} delay';
alignments.s2 = 'T_{2}';
alignments.post_s2 = 'post-T_{2}-delay';
alignments.go = 'go cue';

% iterate through task epochs
task_epochs = fieldnames(rois);
n_epochs = numel(task_epochs);
for ii = 1 : n_epochs
    epoch = task_epochs{ii};
    n_bins.(epoch) = range(rois.(epoch)) * psthbin;
    time.(epoch) = linspace(rois.(epoch)(1),rois.(epoch)(2),n_bins.(epoch));
    psths.(epoch) = nan(n_bins.(epoch),n_neurons,n_i);
end

%% construct T1-aligned, Ii-split psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    
    % trial selection
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
    
    % pre T1 delay-aligned spike rates
    pre_s1_alignment = ...
        repmat(pre_init_padding,n_trials,1);
    pre_s1_alignment_flags = ...
        valid_time >= pre_s1_alignment + rois.pre_s1(1) & ...
        valid_time < pre_s1_alignment + pre_s1_delay(spike_flags);
    pre_s1_chunk_flags = ...
        valid_time >= pre_s1_alignment + rois.pre_s1(1) & ...
        valid_time < pre_s1_alignment + rois.pre_s1(2);
    pre_s1_spkrates = spike_rates;
    pre_s1_spkrates(~pre_s1_alignment_flags') = nan;
    pre_s1_spkrates = reshape(...
        pre_s1_spkrates(pre_s1_chunk_flags'),[n_bins.pre_s1,n_trials])';
    
    % T1-aligned spike rates
    s1_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags);
    s1_alignment_flags = ...
        valid_time >= s1_alignment + rois.s1(1) & ...
        valid_time < s1_alignment + t1(spike_flags);
    s1_chunk_flags = ...
        valid_time >= s1_alignment + rois.s1(1) & ...
        valid_time < s1_alignment + rois.s1(2);
    s1_spkrates = spike_rates;
    s1_spkrates(~s1_alignment_flags') = nan;
    s1_spkrates = reshape(...
        s1_spkrates(s1_chunk_flags'),[n_bins.s1,n_trials])';
    
    % inter T1-T2 delay-aligned spike rates
    inter_s1s2_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags);
    inter_s1s2_alignment_flags = ...
        valid_time >= inter_s1s2_alignment + rois.inter_s1s2(1) & ...
        valid_time < inter_s1s2_alignment + isi;
    inter_s1s2_chunk_flags = ...
        valid_time >= inter_s1s2_alignment + rois.inter_s1s2(1) & ...
        valid_time < inter_s1s2_alignment + rois.inter_s1s2(2);
    inter_s1s2_spkrates = spike_rates;
    inter_s1s2_spkrates(~inter_s1s2_alignment_flags') = nan;
    inter_s1s2_spkrates = reshape(...
        inter_s1s2_spkrates(inter_s1s2_chunk_flags'),[n_bins.inter_s1s2,n_trials])';
    
    % T2-aligned spike rates
    s2_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi;
    s2_alignment_flags = ...
        valid_time >= s2_alignment + rois.s2(1) & ...
        valid_time < s2_alignment + t2(spike_flags);
    s2_chunk_flags = ...
        valid_time >= s2_alignment + rois.s2(1) & ...
        valid_time < s2_alignment + rois.s2(2);
    s2_spkrates = spike_rates;
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = reshape(...
        s2_spkrates(s2_chunk_flags'),[n_bins.s2,n_trials])';
    
    % post T2 delay-aligned spike rates
    post_s2_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi + ...
        t2(spike_flags);
    post_s2_alignment_flags = ...
        valid_time >= post_s2_alignment + rois.post_s2(1) & ...
        valid_time < post_s2_alignment + post_s2_delay;
    post_s2_chunk_flags = ...
        valid_time >= post_s2_alignment + rois.post_s2(1) & ...
        valid_time < post_s2_alignment + rois.post_s2(2);
    post_s2_spkrates = spike_rates;
    post_s2_spkrates(~post_s2_alignment_flags') = nan;
    post_s2_spkrates = reshape(...
        post_s2_spkrates(post_s2_chunk_flags'),[n_bins.post_s2,n_trials])';
    
    % go-aligned spike rates
    go_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi + ...
        t2(spike_flags) + ...
        post_s2_delay;
    go_alignment_flags = ...
        valid_time >= go_alignment + rois.go(1) & ...
        valid_time < go_alignment + rois.go(2);
    go_chunk_flags = go_alignment_flags;
    go_spkrates = spike_rates;
    go_spkrates(~go_alignment_flags') = nan;
    go_spkrates = reshape(...
        go_spkrates(go_chunk_flags'),[n_bins.go,n_trials])';
    
    % compute mean spike density function
    psths.pre_s1(:,nn,i2_mode_idx) = nanmean(pre_s1_spkrates,1);
    psths.s1(:,nn,i2_mode_idx) = nanmean(s1_spkrates,1);
    psths.inter_s1s2(:,nn,i2_mode_idx) = nanmean(inter_s1s2_spkrates,1);
    psths.s2(:,nn,i2_mode_idx) = nanmean(s2_spkrates,1);
    psths.post_s2(:,nn,i2_mode_idx) = nanmean(post_s2_spkrates,1);
    psths.go(:,nn,i2_mode_idx) = nanmean(go_spkrates,1);
end
a=1
%% generate I2-modulated spike rates

% modulation settings
modulation = log(i_set) ./ log(i_set(i2_mode_idx));
scaling = modulation * 1;
gain = modulation * 1;

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'generating fake rates');
    
    % iterate through intensities
    for ii = 1 : n_i
        if ii == i2_mode_idx
            continue;
        end
        
        % apply temporal scaling
        i2_psths(:,nn,ii) = ...
            interp1(roi_time,i2_psths(:,nn,i2_mode_idx),roi_time*scaling(ii));
        i2_psths(roi_time < 0,nn,ii) = i2_psths(roi_time < 0,nn,i2_mode_idx);
        i2_psths(roi_time >= 0,nn,ii) = gain(ii) * i2_psths(roi_time >= 0,nn,ii);
    end
end

%% generate I2-modulated spike trains

% preallocation
if isfield(data,'FakeFR')
    data = rmfield(data,'FakeFR');
end
data.FakeFR = nan(size(data.FR));

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'generating fake spikes');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);

    figure;
    hold on;
    
    % iterate through intensities
    for ii = 1 : n_i
        i2_flags = i2 == i_set(ii);
        i2_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i2_flags;
        spike_trials = find(i2_spike_flags);
        i2_n_trials = numel(spike_trials);

        plot(roi_time,i2_psths(:,nn,ii),...
            'color',i2_clrs(ii,:));
        
        % iterate through trials
        for kk = 1 : i2_n_trials
            trial_idx = spike_trials(kk);
            time_flags = roi_time <= t2(trial_idx);
            lambda = i2_psths(time_flags,nn,i2(trial_idx)==i_set);
            dur = (abs(roi(1)) + t2(trial_idx)) / 1e3;
            [n,ts] = poissonprocess(lambda,dur);
            offset = ...
                pre_init_padding + ...
                pre_s1_delay(trial_idx) + ...
                t1(trial_idx) + ...
                isi + ...
                roi(1);
            data.FakeFR(trial_idx,:) = ...
                histcounts(ts*1e3,[padded_time,n_paddedtimebins]-offset);
        end
    end
end