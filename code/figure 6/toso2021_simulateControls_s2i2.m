%% initialization
if ~exist('data','var')
    toso2021_main;
end

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
roi_window = [0,t_set(end)] + gamma_kernel.paddx;
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
            padded_time >= s2_alignment + roi_window(1) & ...
            padded_time < s2_alignment + roi_window(2);
        s2_chunk_flags = s2_alignment_flags;
        s2_spkrates = spike_rates';
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(...
            s2_spkrates(s2_chunk_flags'),[roi_n_bins,n_trials])';
        
        % compute mean spike density function
        real_psths(:,nn,tt) = nanmean(s2_spkrates,1);
        
%         if ismember(nn,flagged_neurons)
%             figure('position',[1.8,41.8,766.4,740.8]); 
%             subplot(3,1,[1,2]);
%             imagesc(s2_spkrates)
%             subplot(3,1,3);
%             hold on;
%             for ii = 1 : n_epochs
%                 epoch = task_epochs{ii};
%                 plot(roi_time,real_psths(:,nn,tt),...
%                     'linewidth',1.5);
%             end
%             close;
%         end
    end
end

%% construct fake real_psths (for + & - controls)

% iterate through controls
for cc = 1 : n_controls
    control = control_labels{cc};
    
    % control- & feature-specific modulation
    gain = control_modulation.(control).gain;
    offset = control_modulation.(control).offset;
    scaling = control_modulation.(control).scaling;
    
    % iterate through neurons
    for nn = 1 : n_neurons_total
        progressreport(nn,n_neurons_total,...
            sprintf('generating %s control rates',control));
        if all(isnan(real_psths.s2(:,nn)))
            continue;
        end
        
        % iterate through S2 durations
        for tt = 1 : n_t
            
            % iterate through S2 intensities
            for ii = 1 : n_i

                % apply I2 modulation during S2 presentation
                fake_psths.(control)(:,nn,ii) = offset(ii) + gain(ii) * ...
                    interp1(roi_time.s2,...
                    real_psths.s2(:,nn),...
                    roi_time.s2*scaling(ii),...
                    'linear','extrap');
            end
            
            if ismember(nn,flagged_neurons)
                figure;
                hold on;
                for ii = 1 : n_i
                    for ee = 1 : n_epochs
                        epoch = task_epochs{ee};
                        plot(time+offsets,...
                            fake_psths.(control)(:,nn,ii),...
                            'color',i2_clrs(ii,:),...
                            'linewidth',1.5);
                    end
                end
                close;
            end
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
    spike_trials = find(trial_flags);
    n_trials = numel(spike_trials);
    
    %
    if ismember(nn,flagged_neurons)
        figure; hold on;
    end
    
    % iterate through trials
    for kk = 1 : n_trials
        trial_idx = spike_trials(kk);

        %
        lambda_isi = real_psths.isi(:,nn);
        
        % positive control for S2
        lambda_s2 = ...
            fake_psths.positive.s2(roi_time.s2<=t2(trial_idx),nn,i2(trial_idx)==i_set) + ...
            lambda_isi(end) - fake_psths.s2.positive(1,nn,i2(trial_idx)==i_set);
        
        %
        lambda_post_s2 = real_psths.post_s2(:,nn) + ...
            lambda_s2(end) - real_psths.post_s2(1,nn);
        
        lambda = [...
            lambda_pre_s1;...
            lambda_s1;...
            lambda_isi;...
            lambda_s2;...
            lambda_post_s2;...
            lambda_go];
        
        if ismember(nn,flagged_neurons) && t2(trial_idx) == t_set(end)
            
            plot(lambda_s2,...*i2(trial_idx)/80,...
                'color',i2_clrs(i_set==i2(trial_idx),:));
            a=1
        end
        
        % sample spike times
        dur = numel(lambda) * psthbin;
        [n,ts] = poissonprocess(lambda,dur / 1e3);
        offset = pre_init_padding + roi_window.pre_s1(1);
        data.FakeFR(trial_idx,:) = ...
            histcounts(ts*1e3,[gauss_padded_time,n_paddedtimebins*psthbin]-offset);
    end
end