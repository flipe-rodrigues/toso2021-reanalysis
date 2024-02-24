%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% ROI settings

% roi definition
pre_padd = 0;
roi2use = [-pre_padd,t_set(end)];
roi2plot = [-pre_padd,t_set(end)];
roi2plot_padded = round(roi2plot + [-1,1] * .05 * range(roi2plot));
roi2use_n_bins = range(roi2use) / psthbin;
roi2plot_n_bins = range(roi2plot_padded) / psthbin;
roi2use_time = linspace(roi2use(1),roi2use(2),roi2use_n_bins);
roi2plot_time = linspace(roi2plot_padded(1),roi2plot_padded(2),roi2plot_n_bins);
roi2use_flags = ...
    roi2plot_time >= roi2use(1) & ...
    roi2plot_time <= roi2use(2);

%% subject selection
subject_flags = ismember(subjects,subject_set);

%% exponential decay function
taufun = @(tau,x) exp(-x * tau);

%% construct s2-aligned psths

% preallocation
R = nan(roi2plot_n_bins,n_neurons,n_contrasts);
taus = nan(n_neurons,n_contrasts);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'computing average activity');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    ref_trial_flags = ...
        valid_flags & ...
        subject_flags & ...
        neuron_flags;
    if sum(ref_trial_flags) == 0
        continue;
    end
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        trial_flags = ...
            valid_flags & ...
            subject_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(trial_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_rates = data.SDF(trial_flags,:);
        n_trials = size(spike_rates,1);

        % S2-aligned spike rates
        alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(trial_flags) + ...
            t1(trial_flags) + ...
            isi;
        alignment_flags = ...
            padded_time >= alignment_onset + roi2plot_padded(1) & ...
            padded_time < alignment_onset + t2(trial_flags);
        chunk_flags = ...
            padded_time >= alignment_onset + roi2plot_padded(1) & ...
            padded_time < alignment_onset + roi2plot_padded(2);
        spkrates = spike_rates';
        spkrates(~alignment_flags') = nan;
        spkrates = reshape(...
            spkrates(chunk_flags'),[roi2plot_n_bins,n_trials])';

        % compute autocorrelogram
        spkrate = nanmean(spkrates,1);
        nan_flags = isnan(spkrate);
        spkrate(nan_flags) = 0;
        [r,lags] = xcorr(spkrate,spkrate);
        lag_flags = lags >= 0;
        [fitobj,gof] = fit(lags(lag_flags)'/1e3,r(lag_flags)'/max(r),taufun,...
            'lower',0,...
            'upper',inf,...
            'start',1);
        taus(nn,ii) = coeffvalues(fitobj);
        R(:,nn,ii) = r(lag_flags)'/max(r);
    end
end

%%
figure;
hold on;

[~,edges] = histcounts(taus);
edges = linspace(1,3,50);

for ii = n_contrasts : -1 : 1
    histogram(taus(:,ii),edges,...
        'normalization','pdf',...
        'facecolor',contrast_clrs(ii,:));
end

figure;
hold on;
for ii = n_contrasts : -1 : 1
    plot(lags(lag_flags),nanmean(R(:,:,ii),2),...
        'color',contrast_clrs(ii,:),...
        'linewidth',1.5);
end