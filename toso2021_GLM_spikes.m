%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% GLM settings
distro = 'normal';
glm_roi = [0,150] - t2_set(1);
spkintegration_window = diff(glm_roi);

%% construct response

% preallocation
spkcounts = nan(n_total_trials,1);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'fetching spike counts for GLM');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % flag trials for the current condition
    s2_spike_flags = ...
        valid_flags & ...
        neuron_flags;
    flagged_trials = find(s2_spike_flags);
    if sum(s2_spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    s2_spike_counts = data.FR(s2_spike_flags,:);
    s2_spike_rates = ...
        conv2(1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
    n_trials = size(s2_spike_counts,1);
    
    % T2-offset-aligned spike rates
    s2_alignment = ...
        pre_init_padding + ...
        pre_t1_delay(s2_spike_flags) + ...
        t1(s2_spike_flags) + ...
        isi + ...
        t2(s2_spike_flags);
    s2_alignment_flags = ...
        valid_time >= s2_alignment + glm_roi(1) & ...
        valid_time < s2_alignment + glm_roi(2);
    s2_chunk_flags = s2_alignment_flags;
    s2_spkrates = s2_spike_rates;
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = ...
        reshape(s2_spkrates(s2_chunk_flags'),[spkintegration_window,n_trials])';

    % store average spike rates
    spkcounts(s2_spike_flags) = nanmean(s2_spkrates,2);
end

%% spike count GLM

% design matrix
X = [s1,s2,d1,d2];
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
        'predictorvars',{s1_lbl,s2_lbl,d1_lbl,d2_lbl},...
        'distribution',distro,...
        'intercept',true);
    betas(nn,:) = mdl.Coefficients.Estimate;
    pvals(nn,:) = mdl.Coefficients.pValue;
end

%% rectification
rectified_betas = (betas);

%% significance
significant_betas = rectified_betas;
significant_betas(pvals>=.05) = nan;
insignificant_betas = rectified_betas;
insignificant_betas(pvals<.05) = nan;

%% plot GLM coefficients
fig = figure(figopt,...
    'name','choice_GLM');
axes(axesopt.default,...
    'xlim',[0,n_coefficients+1],...
    'ylim',[-3,3]+[-1,1]*.05*6,...
    'xtick',1:n_coefficients,...
    'xticklabelrotation',0,...
    'xticklabel',beta_labels,...
    'clipping','off',...
    'plotboxaspectratio',[1,1,1]);
title(sprintf('%s>%s~%s(\\phi(\\betaX))',s2_lbl,s1_lbl,capitalize(distro)));
xlabel('Regressor');
ylabel('|Weight|');

% iterate through coefficients
for bb = 1 : n_coefficients
    offset = ((1:n_neurons) - (n_neurons + 1) / 2) * .0001 * range(xlim);
   
    % plot subject coefficients
    grapeplot(bb+offset,insignificant_betas(:,bb),...
        'marker','o',...
        'markersize',5,...
        'markerfacecolor',subject_clr,...
        'markeredgecolor','w',...
        'linewidth',1.5);
    grapeplot(bb+offset,significant_betas(:,bb),...
        'marker','o',...
        'markersize',5,...
        'markerfacecolor','k',...
        'markeredgecolor','w',...
        'linewidth',1.5);
end

% plot animal-pool coefficients
p = stem(1:n_coefficients,nanmean(rectified_betas,1),...
    'color','k',...
    'marker','o',...
    'markersize',12,...
    'markerfacecolor','k',...
    'markeredgecolor','w',...
    'linewidth',1.5);
p.BaseLine.LineWidth = p.LineWidth;

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end