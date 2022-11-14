%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct T1-aligned, Ii-split psths
roi = [-500,t_set(end)];
roi_n_bins = range(roi) * psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

% preallocation
i1_psths = nan(roi2plot_n_bins,n_neurons,n_i);
i2_psths = nan(roi2plot_n_bins,n_neurons,n_i);

% smoothing settings
gauss_kernel = gausskernel('mu',0,'sig',50,'binwidth',psthbin);
    
% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    i1_flags = i1 == i_set(i2_mode_idx);
    i2_flags = i2 == i_set(i2_mode_idx);
    i1_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i1_flags;
    i2_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i2_flags;
    if sum(i1_spike_flags) == 0 || ...
            sum(i2_spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    i2_spike_counts = data.FR(i2_spike_flags,:);
    i2_spike_rates = conv2(...
        1,gauss_kernel.pdf,i2_spike_counts,'valid')' / psthbin * 1e3;
    i2_n_trials = size(i2_spike_counts,1);
    
    % S2-aligned spike rates
    alignment_onset = ...
        pre_init_padding + ...
        pre_t1_delay(i2_spike_flags) + ...
        t1(i2_spike_flags) + ...
        isi;
    alignment_flags = ...
        valid_time >= alignment_onset + roi2plot(1) & ...
        valid_time < alignment_onset + t2(i2_spike_flags);
    chunk_flags = ...
        valid_time >= alignment_onset + roi2plot(1) & ...
        valid_time < alignment_onset + roi2plot(2);
    spkrates = i2_spike_rates;
    spkrates(~alignment_flags') = nan;
    spkrates = reshape(...
        spkrates(chunk_flags'),[roi2plot_n_bins,i2_n_trials])';
    
    % compute mean spike density function
    i2_psths(:,nn,i2_mode_idx) = nanmean(spkrates,1);
    
    % fetch spike counts & compute spike rates
    i2_spike_counts = data.FR(i2_spike_flags,:);
    i2_spike_rates = conv2(...
        1,gauss_kernel.pdf,i2_spike_counts,'valid')' / psthbin * 1e3;
    i2_n_trials = size(i2_spike_counts,1);
    
    % S2-aligned spike rates
    alignment_onset = ...
        pre_init_padding + ...
        pre_t1_delay(i2_spike_flags) + ...
        t1(i2_spike_flags) + ...
        isi;
    alignment_flags = ...
        valid_time >= alignment_onset + roi2plot(1) & ...
        valid_time < alignment_onset + t2(i2_spike_flags);
    chunk_flags = ...
        valid_time >= alignment_onset + roi2plot(1) & ...
        valid_time < alignment_onset + roi2plot(2);
    spkrates = i2_spike_rates;
    spkrates(~alignment_flags') = nan;
    spkrates = reshape(...
        spkrates(chunk_flags'),[roi2plot_n_bins,i2_n_trials])';
    
    % compute mean spike density function
    i2_psths(:,nn,i2_mode_idx) = nanmean(spkrates,1);
end

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