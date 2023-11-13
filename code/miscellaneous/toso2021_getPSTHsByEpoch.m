%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% epoch settings

% preallocation
epochs = struct();
psths = struct();
zpsths = struct();

% epoch alignments
epochs.label.pre_s1 = 'pre-T_{1}-delay';
epochs.label.s1 = 'T_{1}';
epochs.label.isi = 'ISI';
epochs.label.s2 = 'T_{2}';
epochs.label.post_s2 = 'post-T_{2}-delay';
epochs.label.go = 'go cue';

% epoch rois
epochs.roi.pre_s1 = [-500,unique(pre_t1_delay(valid_flags))];
epochs.roi.s1 = [0,t_set(end)];
epochs.roi.isi = [0,isi];
epochs.roi.s2 = [0,t_set(end)];
epochs.roi.post_s2 = [0,1] * post_t2_delay;
epochs.roi.go = [0,450];

% epoch contrasts
epochs.contrast.pre_s1 = 'prev_t2';
epochs.contrast.s1 = 'i1';
epochs.contrast.isi = 't1';
epochs.contrast.s2 = 'i2';
epochs.contrast.post_s2 = 'choice';
epochs.contrast.go = 'choice';

% iterate through epochs
epoch_labels = fieldnames(epochs.label);
n_epochs = numel(epoch_labels);
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    epoch_contrast_str = epochs.contrast.(epoch);
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    epochs.n_bins.(epoch) = range(epochs.roi.(epoch)) * psthbin;
    epochs.time.(epoch) = linspace(epochs.roi.(epoch)(1),epochs.roi.(epoch)(2),epochs.n_bins.(epoch));
    psths.(epoch) = nan(epochs.n_bins.(epoch),n_neurons,epoch_n_contrasts);
end

%% for cross-epoch intersectional selection
surviving_neurons = true(n_neurons,1);

%% construct pre-S1-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing pre-s1 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.pre_s1;
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str(1:end-1),'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
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
            valid_time >= pre_s1_alignment + epochs.roi.pre_s1(1) & ...
            valid_time < pre_s1_alignment + pre_t1_delay(spike_flags);
        pre_s1_chunk_flags = ...
            valid_time >= pre_s1_alignment + epochs.roi.pre_s1(1) & ...
            valid_time < pre_s1_alignment + epochs.roi.pre_s1(2);
        pre_s1_spkrates = spike_rates;
        pre_s1_spkrates(~pre_s1_alignment_flags') = nan;
        pre_s1_spkrates = reshape(...
            pre_s1_spkrates(pre_s1_chunk_flags'),[epochs.n_bins.pre_s1,n_trials])';
        
        % compute mean spike density function
        psths.pre_s1(:,nn,ii) = nanmean(pre_s1_spkrates,1);
    end
end

%% construct S1-onset-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing s1 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.s1;
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str(1:end-1),'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % T1-aligned spike rates
        s1_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags);
        s1_alignment_flags = ...
            valid_time >= s1_alignment + epochs.roi.s1(1) & ...
            valid_time < s1_alignment + t1(spike_flags);
        s1_chunk_flags = ...
            valid_time >= s1_alignment + epochs.roi.s1(1) & ...
            valid_time < s1_alignment + epochs.roi.s1(2);
        s1_spkrates = spike_rates;
        s1_spkrates(~s1_alignment_flags') = nan;
        s1_spkrates = reshape(...
            s1_spkrates(s1_chunk_flags'),[epochs.n_bins.s1,n_trials])';
        
        % compute mean spike density function
        psths.s1(:,nn,ii) = nanmean(s1_spkrates,1);
    end
end

%% construct ISI-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing isi psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.isi;
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str(1:end-1),'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % ISI-aligned spike rates
        isi_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags);
        isi_alignment_flags = ...
            valid_time >= isi_alignment + epochs.roi.isi(1) & ...
            valid_time < isi_alignment + isi;
        isi_chunk_flags = ...
            valid_time >= isi_alignment + epochs.roi.isi(1) & ...
            valid_time < isi_alignment + epochs.roi.isi(2);
        isi_spkrates = spike_rates;
        isi_spkrates(~isi_alignment_flags') = nan;
        isi_spkrates = reshape(...
            isi_spkrates(isi_chunk_flags'),[epochs.n_bins.isi,n_trials])';
        
        % compute mean spike density function
        psths.isi(:,nn,ii) = nanmean(isi_spkrates,1);
    end
end

%% construct S2-onset-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing s2 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.s2;
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str(1:end-1),'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % T2-aligned spike rates
        s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi;
        s2_alignment_flags = ...
            valid_time >= s2_alignment + epochs.roi.s2(1) & ...
            valid_time < s2_alignment + t2(spike_flags);
        s2_chunk_flags = ...
            valid_time >= s2_alignment + epochs.roi.s2(1) & ...
            valid_time < s2_alignment + epochs.roi.s2(2);
        s2_spkrates = spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(...
            s2_spkrates(s2_chunk_flags'),[epochs.n_bins.s2,n_trials])';
        
        % compute mean spike density function
        psths.s2(:,nn,ii) = nanmean(s2_spkrates,1);
    end
end

%% construct post-S2-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing post-s2 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.post_s2;
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str(1:end-1),'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % post T2 delay-aligned spike rates
        post_s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi + ...
            t2(spike_flags);
        post_s2_alignment_flags = ...
            valid_time >= post_s2_alignment + epochs.roi.post_s2(1) & ...
            valid_time < post_s2_alignment + post_t2_delay;
        post_s2_chunk_flags = ...
            valid_time >= post_s2_alignment + epochs.roi.post_s2(1) & ...
            valid_time < post_s2_alignment + epochs.roi.post_s2(2);
        post_s2_spkrates = spike_rates;
        post_s2_spkrates(~post_s2_alignment_flags') = nan;
        post_s2_spkrates = reshape(...
            post_s2_spkrates(post_s2_chunk_flags'),[epochs.n_bins.post_s2,n_trials])';
        
        % compute mean spike density function
        psths.post_s2(:,nn,ii) = nanmean(post_s2_spkrates,1);
    end
end

%% construct go-cue-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing go psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.go;
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str(1:end-1),'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % go-aligned spike rates
        go_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi + ...
            t2(spike_flags) + ...
            post_t2_delay;
        go_alignment_flags = ...
            valid_time >= go_alignment + epochs.roi.go(1) & ...
            valid_time < go_alignment + epochs.roi.go(2);
        go_chunk_flags = go_alignment_flags;
        go_spkrates = spike_rates;
        go_spkrates(~go_alignment_flags') = nan;
        go_spkrates = reshape(...
            go_spkrates(go_chunk_flags'),[epochs.n_bins.go,n_trials])';
        
        % compute mean spike density function
        psths.go(:,nn,ii) = nanmean(go_spkrates,1);
    end
end

%% intersectional neuron selection
psths.pre_s1 = psths.pre_s1(:,surviving_neurons,:);
psths.s1 = psths.s1(:,surviving_neurons,:);
psths.isi = psths.isi(:,surviving_neurons,:);
psths.s2 = psths.s2(:,surviving_neurons,:);
psths.post_s2 = psths.post_s2(:,surviving_neurons,:);
psths.go = psths.go(:,surviving_neurons,:);
n_neurons = sum(surviving_neurons);

%% normalization
total_n_bins = 0;

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    total_n_bins = total_n_bins + epochs.n_bins.(epoch);
end

% preallocation
psth_concat = nan(total_n_bins,n_neurons);

% iterate through task epochs
last_idx = 0;
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    
    % iterate through contrasts
    for jj = 1 : epoch_n_contrasts
        idcs = (1 : epochs.n_bins.(epoch)) + last_idx;
        psth_concat(idcs,:) = psths.(epoch)(:,:,jj);
        last_idx = idcs(end);
    end
end

% z-scoring
mus = nanmean(psth_concat,1);
sigs = nanstd(psth_concat,0,1);
zpsth_concat = (psth_concat - mus) ./ sigs;

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    
    % z-scoring
    zpsths.(epoch) = (psths.(epoch) - mus) ./ sigs;
end