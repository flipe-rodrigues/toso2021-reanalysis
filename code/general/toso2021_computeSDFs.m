%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% compute spike density functions

% check if SDFs have already been computed
if ~isfield(data,'SDF')
    
    % preallocation
    data.SDF = nan(size(data.FR));
    
    % iterate through neurons
    for nn = 1 : n_neurons_total
        progressreport(nn,n_neurons_total,...
            'computing spike density functions');
        trial_flags = data.NeuronNumb == neuron_idcs(nn);
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(trial_flags,:);
        spike_rates = conv2(...
            1,gamma_kernel.pdf,spike_counts,'valid') / psthbin * 1e3;

        % store spike density function
        data.SDF(trial_flags,validtime_flags) = spike_rates;
    end
end