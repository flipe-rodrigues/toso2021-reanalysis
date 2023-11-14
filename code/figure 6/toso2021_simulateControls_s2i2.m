%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% modulation settings
modulation = log(i_set') ./ log(i_set(i2_mode_idx));

%% control settings
control_modulation = struct(...
    'positive',struct(...
        'gain',modulation.^1,...
        'offset',modulation.*0,...
        'scaling',modulation.^0),...
    'negative',struct(...
        'gain',modulation.^0,...
        'offset',modulation.*0,...
        'scaling',modulation.^0));
control_labels = fieldnames(control_modulation);
n_controls = numel(control_labels);

% iterate through controls
for cc = 1 : n_controls
    control = control_labels{cc};
    
    % display control settings
    fprintf('%s control modulation\n\t\t  I2 = {%i, %i, %i} %s:\n',...
        upper(control),i_set,i2_units);
    disp(control_modulation.positive);
end

%% ROI settings
roi_padding = gamma_kernel.paddx * 1;
roi_window = [0,t_set(end)] + roi_padding;
roi_n_bins = range(roi_window) / psthbin;
roi_time = linspace(roi_window(1),roi_window(2),roi_n_bins);

%% construct real real_psths from middle-I2 trials

% preallocation
real_psths = nan(roi_n_bins,n_neurons_total,n_t);
fake_psths = struct(...
    control_labels{1},nan(roi_n_bins,n_neurons_total,n_t,n_i),...
    control_labels{2},nan(roi_n_bins,n_neurons_total,n_t,n_i));

% iterate through neurons
for nn = 1 : n_neurons_total
    progressreport(nn,n_neurons_total,'parsing neural data');
    neuron_flags = data.NeuronNumb == nn;
    
    % iterate through S2 durations
    for tt = 1 : n_t
        t2_flags = t2 == t_set(tt);
        
        % I2-based trial selection
        i2_flags = i2 == i_set(i2_mode_idx);
        
        % trial selection
        trial_flags = ...
            valid_flags & ...
            neuron_flags & ...
            t2_flags & ...
            i2_flags;
        n_trials = sum(trial_flags);
        if n_trials == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_rates = data.SDF(trial_flags,:);
        
        % S2-aligned spike rates
        s2_alignment = ...
            pre_init_padding + ...
            pre_s1_delay(trial_flags) + ...
            t1(trial_flags) + ...
            isi;
        s2_alignment_flags = ...
            padded_time >= s2_alignment + roi_padding(1) & ...
            padded_time < s2_alignment + t2(trial_flags) + roi_padding(2);
        s2_chunk_flags = ...
            padded_time >= s2_alignment + roi_window(1) & ...
            padded_time < s2_alignment + roi_window(2);
        s2_spkrates = spike_rates';
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(...
            s2_spkrates(s2_chunk_flags'),[roi_n_bins,n_trials])';
        
        % reconvolve with a flipped kernel
%         s2_spkrates = conv2(1,fliplr(gamma_kernel.pdf),s2_spkrates,'valid');
        
        % compute mean spike density function
        real_psths(:,nn,tt) = nanmean(s2_spkrates,1);
        
%         if ismember(nn,flagged_neurons)
%             figure('position',[1.8,41.8,766.4,740.8]); 
%             subplot(3,1,[1,2]);
%             imagesc(roi_window,[],s2_spkrates)
%             subplot(3,1,3);
%             plot(roi_time,real_psths(:,nn,tt),...
%                 'linewidth',1.5);
%             xlim(roi_window);
%             close;
%         end
    end
end

%% construct fake real_psths (for + & - controls)

% iterate through controls
for cc = 1 : n_controls
    control = control_labels{cc};
    
    % control- & feature-specific modulation
    nominal_gain = control_modulation.(control).gain;
    
    % iterate through neurons
    for nn = 1 : n_neurons_total
        progressreport(nn,n_neurons_total,...
            sprintf('generating %s control rates',control));
        
        % iterate through S2 durations
        for tt = 1 : n_t
            s2_time_flags = ...
                roi_time >= 0 & ...
                roi_time < t_set(tt);
            post_s2_time_flags = ...
                roi_time >= t_set(tt);

            % iterate through S2 intensities
            for ii = 1 : n_i

                % compute gain through time
                gain_through_time = ones(roi_n_bins,1);
                gain_through_time(s2_time_flags) = nominal_gain(ii);
                gain_through_time(post_s2_time_flags) = ...
                    linspace(nominal_gain(ii),1,sum(post_s2_time_flags));
                
                % apply I2 modulation during S2 presentation
                fake_psths.(control)(:,nn,tt,ii) = ...
                    gain_through_time .* real_psths(:,nn,tt);
            end
            
%             if ismember(nn,flagged_neurons)
%                 figure;
%                 hold on;
%                 for ii = 1 : n_i
%                     plot(roi_time,fake_psths.(control)(:,nn,tt,ii),...
%                         'color',i2_clrs(ii,:),...
%                         'linewidth',1.5);
%                 end
%                 plot(roi_time,gain_through_time*5,'k');
%                 close;
%             end
        end
    end
end

%% sample fake spike times (for + & - controls)

% preallocation
if ~isfield(data,'FakeFR')
    data.FakeFR = nan(size(data.FR));
end

% iterate through neurons
for nn = 1 : n_neurons_total
    progressreport(nn,n_neurons_total,'generating fake spikes');
    neuron_flags = data.NeuronNumb == nn;
    
    % trial selection
    trial_flags = ...
        valid_flags & ...
        neuron_flags;
    trial_idcs = find(trial_flags);
    n_trials = numel(trial_idcs);
    
    %
%     if ismember(nn,flagged_neurons)
%         figure; hold on;
%     end
    
    % iterate through trials
    for kk = 1 : n_trials
        trial_idx = trial_idcs(kk);

        % preallocation
        lambda = nan(n_paddedtimebins,1);
        
        % S1-alignment
        s1_trial_alignment = ...
            pre_init_padding + ...
            pre_s1_delay(trial_idx);
        s1_trial_alignment_flags = ...
            padded_time >= s1_trial_alignment + roi_window(1) & ...
            padded_time < s1_trial_alignment + roi_window(2);
        
        % S2-alignment
        s2_trial_alignment = ...
            pre_init_padding + ...
            pre_s1_delay(trial_idx) + ...
            t1(trial_idx) + ...
            isi;
        s2_trial_alignment_flags = ...
            padded_time >= s2_trial_alignment + roi_window(1) & ...
            padded_time < s2_trial_alignment + roi_window(2);

        % negative control for S1
        lambda(s1_trial_alignment_flags) = fake_psths.negative(:,nn,...
            t1(trial_idx)==t_set,i2(trial_idx)==i_set);
        
        % positive control for S2
        lambda(s2_trial_alignment_flags) = fake_psths.positive(:,nn,...
            t2(trial_idx)==t_set,i2(trial_idx)==i_set);

%         if ismember(nn,flagged_neurons)
%             plot(padded_time-padded_time(s1_trial_alignment),lambda,...
%                 'color',i2_clrs(i_set==i2(trial_idx),:),...
%                 'linewidth',1.5);
%             title(sprintf('T1 = %i; T2 = %i',t1(trial_idx),t2(trial_idx)));
%             [t1(trial_idx),t2(trial_idx)]
%         end
        
        % sample spike times
        dur = numel(lambda) * psthbin;
        [n,ts] = poissonprocess(lambda,dur / 1e3);
        fake_spike_counts = histcounts(...
            ts*1e3,[padded_time,n_paddedtimebins*psthbin]);
        data.FakeFR(trial_idx,:) = fake_spike_counts;
    end
end