%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% GLM settings
distro = 'poisson';
glm_roi = [0,t_set(t2_mode_idx)];
glm_win = diff(glm_roi);

%% construct response

% preallocation
spkcounts = nan(n_total_trials,1);

% I1 clamping
i1_flags = i1 == i_set(i1_mode_idx);

% T2 selection
t2_flags = t2 >= glm_roi(2);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'fetching spike counts for GLM');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % flag trials for the current condition
    spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        ...i1_flags & ...
        t2_flags;
    flagged_trials = find(spike_flags);
    if sum(spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    spike_counts = data.FR(spike_flags,:)';
    spike_rates = ...
        conv2(kernel.pdf,1,spike_counts,'valid') / psthbin * 1e3;
    n_trials = sum(spike_flags);
    
    % T2-offset-aligned spike rates
    alignment_onset = ...
        pre_init_padding + ...
        pre_t1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi;
    alignment_flags = ...
        padded_time >= alignment_onset + glm_roi(1) & ...
        padded_time < alignment_onset + glm_roi(2);
    chunk_flags = alignment_flags;
    spkrates = spike_counts;
    spkrates(~alignment_flags') = nan;
    spkrates = reshape(spkrates(chunk_flags'),[glm_win,n_trials])';
    
    % store average spike rates
    spkcounts(spike_flags) = nansum(spkrates,2);
end

%% spike count GLM

% design matrix
X = [d1,d2];
n_regressors = size(X,2);
n_coefficients = n_regressors + 1;

% feature normalization
Z = (X - nanmean(X)) ./ nanstd(X);

% preallocation
betas = nan(n_neurons,n_coefficients);
pvals = nan(n_neurons,n_coefficients);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'fitting neuron-wise GLM');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    trial_flags = ...
        valid_flags & ...
        neuron_flags;
    
    % fit GLM to each subject
    mdl = fitglm(Z(trial_flags,:),spkcounts(trial_flags),'linear',...
        'predictorvars',{d1_lbl,d2_lbl},...
        'distribution',distro,...
        'intercept',true);
    betas(nn,:) = mdl.Coefficients.Estimate;
    pvals(nn,:) = mdl.Coefficients.pValue;
end

%% print percentage of significantly modulated neurons
alpha = .05;
significance_mask = squeeze(pvals(:,2:end)) < alpha;
n_significant = nansum(significance_mask);

% iterate through coefficients
for bb = 1 : n_regressors
    fprintf('%s-modulated: %i/%i (%.1f%%)\n',...
        mdl.Coefficients.Properties.RowNames{bb+1},...
        n_significant(bb),n_neurons,...
        n_significant(bb)/n_neurons*100);
end

%% plot spike count distribution

% bin settings
n_bins = 21;
binspan = [0,20];
binedges = linspace(binspan(1),binspan(2),n_bins);
egneuron_flags = ismember(data.NeuronNumb,...
    [68,72,215,391,393,428,459,470,526]);
bincounts = histcounts(spkcounts(valid_flags),... & egneuron_flags),...
    'binedges',binedges);
bincounts = bincounts / nansum(bincounts);

% figure initialization
if glm_roi(2) ~= 100
    fig = figure(figopt,...
        'name','GLM_spikeCountDistro');
end
set(gca,...
    axesopt.default,...
    'xlim',binspan+[-1,1]*.05*range(binspan),...
    'ylimspec','tight',...
    'ytick',0,...
    'ycolor','k',...
    'clipping','off',...
    'plotboxaspectratio',[2.25,1,1]);
xlabel('Spike count');
ylabel('PDF');

% plot spike count distribution
if glm_roi(2) == 100
    histogram(...
        'binedges',binedges,...
        'bincounts',bincounts,...
        'facealpha',.75,...
        'edgecolor','k',...
        'facecolor','w',...
        'linewidth',1.5);
else
    histogram(...
        'binedges',binedges,...
        'bincounts',bincounts,...
        'facealpha',1,...
        'edgecolor','w',...
        'facecolor','k',...
        'linewidth',1.5);
end

% update y-axis
ylim(ylim+[0,1]*.05*range(ylim));

% legend
% leg = legend({'334 ms','100 ms'},...
%     'position',[0.085,0.66,.27,.2],...
%     'box','off');
% title(leg,'Spike integration window:');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end