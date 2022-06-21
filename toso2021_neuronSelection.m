%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% manual curation

% selected for being unstable (assessed by looking at spike rasters)
if strcmpi(task_str,'duration')
    neurons2exclude = [...
        1,19,29,51,60,63,71,86,92,...
        153,239,326,453,...
        557,584,585,599];
elseif strcmpi(task_str,'intensity')
   neurons2exclude = []; 
end

%% construct S2-aligned, Ti- & Ii-split psths

% time settings
roi = [0,t_set(end)];
roi_n_bins = range(roi) * psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

% preallocation
mean_frs = nan(numel(neuron_idcs),n_t,n_i);
trial_type_counts = nan(numel(neuron_idcs),n_t,n_i);
stability_coeffs = nan(numel(neuron_idcs),1);

% iterate through neurons
for nn = neuron_idcs'
    progressreport(nn,numel(neuron_idcs),'computing selection criteria');
    neuron_flags = data.NeuronNumb == neuron_idcs(nn);
    
    % iterate through durations
    for tt = 1 : n_t
        t2_flags = t2 == t_set(tt);
        
        % iterate through intensities
        for ii = 1 : n_i
            i1_flags = i1 == i_set(ii);
            i2_flags = i2 == i_set(ii);
            s2i1_spike_flags = ...
                valid_flags & ...
                neuron_flags & ...
                t2_flags & ...
                i1_flags;
            s2i2_spike_flags = ...
                valid_flags & ...
                neuron_flags & ...
                t2_flags & ...
                i2_flags;
            if (sum(s2i1_spike_flags) == 0) || (sum(s2i2_spike_flags) == 0) 
                continue;
            end
            
            % fetch spike counts & compute spike rates
            s2_spike_counts = data.FR(s2i2_spike_flags,:);
            s2_spike_rates = conv2(...
                1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
            s2_n_trials = size(s2_spike_counts,1);
            
            % T2-aligned spike rates
            s2_alignment_offset = ...
                pre_init_padding + ...
                pre_t1_delay(s2i2_spike_flags) + ...
                t1(s2i2_spike_flags) + ...
                inter_t1t2_delay;
            s2_alignment_flags = ...
                valid_time >= s2_alignment_offset + roi(1) & ...
                valid_time < s2_alignment_offset + t2(s2i2_spike_flags);
            s2_chunk_flags = ...
                valid_time >= s2_alignment_offset + roi(1) & ...
                valid_time < s2_alignment_offset + roi(2);
            s2_spkrates = s2_spike_rates;
            s2_spkrates(~s2_alignment_flags') = nan;
            s2_spkrates = reshape(...
                s2_spkrates(s2_chunk_flags'),[roi_n_bins,s2_n_trials])';
            
            % neuron selection criteria
            mean_frs(nn,tt,ii) = nanmean(s2_spkrates,[1,2]);
            trial_type_counts(nn,tt,ii) = s2_n_trials;
        end  
    end
    
    % for computing stability
    i2_flags = i2 == i_set(i2_mode_idx);
    s2i2_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i2_flags;
    if sum(s2i2_spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    s2_spike_counts = data.FR(s2i2_spike_flags,:);
    s2_spike_rates = conv2(...
        1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
    s2_n_trials = size(s2_spike_counts,1);
    
    % T2-aligned spike rates
    s2_alignment_offset = ...
        pre_init_padding + ...
        pre_t1_delay(s2i2_spike_flags) + ...
        t1(s2i2_spike_flags) + ...
        inter_t1t2_delay;
    s2_alignment_flags = ...
        valid_time >= s2_alignment_offset + roi(1) & ...
        valid_time < s2_alignment_offset + t2(s2i2_spike_flags);
    s2_chunk_flags = ...
        valid_time >= s2_alignment_offset + roi(1) & ...
        valid_time < s2_alignment_offset + t_set(t2_mode_idx);
    s2_spkrates = s2_spike_rates;
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = reshape(...
        s2_spkrates(s2_chunk_flags'),[t_set(t2_mode_idx)*psthbin,s2_n_trials])';
    
    % compute stability coefficient
    first_third_idcs = 1 : round(s2_n_trials * 1 / 3);
    last_third_idcs = round(s2_n_trials * 2 / 3) : s2_n_trials;
    first_third_mu = nanmean(s2_spkrates(first_third_idcs,:),1)';
    last_third_mu = nanmean(s2_spkrates(last_third_idcs,:),1)';
    corr_coeffs = corrcoef(...
        first_third_mu(~isnan(first_third_mu)&~isnan(last_third_mu)),...
        last_third_mu(~isnan(first_third_mu)&~isnan(last_third_mu)));
    stability_coeffs(nn) = abs(corr_coeffs(1,2));
end
trial_type_counts(isnan(trial_type_counts)) = 0;

%% neuron selection
mean_fr_flags = ...
    mean(mean_frs(:,:,i2_mode_idx),2) >= mean_fr_cutoff;
trial_count_flags = ...
    all(trial_type_counts >= trial_count_cutoff,[2,3]);
stability_flags = ...
    ~ismember(neuron_idcs,neurons2exclude);
neurons2keep_flags = ...
    mean_fr_flags & ...
    trial_count_flags & ...
    stability_flags;
flagged_neurons = neuron_idcs(neurons2keep_flags);
n_neurons = sum(neurons2keep_flags);

% display selection outcome
fprintf('\nSELECTION CRITERIA:\n');
fprintf('- minimum mean firing rate in middle I2 trials: %i Hz\n',...
    mean_fr_cutoff);
fprintf('- minimum trial count on all T2-I2 combinations: %i\n',...
    trial_count_cutoff);
fprintf('- stability assessed by visual inspection\n');
fprintf('%i/%i neurons passed.\n\n',n_neurons,n_neurons_total);