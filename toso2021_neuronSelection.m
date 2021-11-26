%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct S2-aligned, Ti- & Ii-split psths

% time settings
roi = [0,t_set(end)];
roi_n_bins = range(roi) * psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

% preallocation
mean_frs = nan(numel(neuron_idcs),n_t,n_i);
trial_type_numbers = nan(numel(neuron_idcs),n_t,n_i);

% iterate through neurons
for nn = neuron_idcs'
    progressreport(nn,numel(neuron_idcs),'parsing neural data');
    neuron_flags = data.NeuronNumb == neuron_idcs(nn);
    
    % iterate through stimuli
    for tt = 1 : n_t
        t2_flags = t2 == t_set(tt);
        
        % iterate through intensities
        for ii = 1 : n_i
            i2_flags = i2 == i_set(ii);
            s2_spike_flags = ...
                valid_flags & ...
                neuron_flags & ...
                t2_flags & ...
                i2_flags;
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
                valid_time >= s2_alignment_offset + roi(1) & ...
                valid_time < s2_alignment_offset + t2(s2_spike_flags);
            s2_chunk_flags = ...
                valid_time >= s2_alignment_offset + roi(1) & ...
                valid_time < s2_alignment_offset + roi(2);
            s2_spkrates = s2_spike_rates;
            s2_spkrates(~s2_alignment_flags') = nan;
            s2_spkrates = reshape(...
                s2_spkrates(s2_chunk_flags'),[roi_n_bins,s2_n_trials])';
            
            % neuron selection criteria
            mean_frs(nn,tt,ii) = nanmean(s2_spkrates,[1,2]);
            trial_type_numbers(nn,tt,ii) = s2_n_trials;
        end
    end
end

%% neuron selection
n_neurons_total = numel(neuron_idcs);
mean_fr_flags = ...
    all(mean_frs(:,:,i2_mode_idx) >= mean_fr_cutoff,[2,3]);
trial_number_flags = ...
    all(trial_type_numbers >= n_trial_cutoff,[2,3]);
neurons2keep_flags = ...
    mean_fr_flags & ...
    trial_number_flags;
flagged_neurons = neuron_idcs(neurons2keep_flags);
n_neurons = sum(neurons2keep_flags);

% display selection outcome
fprintf('\nSELECTION CRITERIA:\n');
fprintf('- minimum mean firing rate in middle I2 trials: %i Hz\n',...
    mean_fr_cutoff);
fprintf('- minimum trial count on all T2-I2 combinations: %i\n',...
    n_trial_cutoff);
fprintf('%i/%i neurons passed.\n\n',n_neurons,n_neurons_total);