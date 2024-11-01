%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% training & test set conditions

% training set conditions
if strcmpi(contrast_str,'t1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],t_set([1,end]),...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],t_set([1,end]),...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
end

% test set conditions
if strcmpi(contrast_str,'t1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set,[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
end

% print training & test conditions
fprintf('\nTRAINING CONDITIONS:\n');
conditions.train.values
fprintf('\nTEST CONDITIONS:\n');
conditions.test.values

%% run settings
n_runs = 15;

%% condition weights

% iterate through training conditions
for kk = 1 : conditions.train.n
    
    % flag trials for the current condition
    feature_flags = false(n_total_trials,conditions.train.features.n);
    for ff = 1 : conditions.train.features.n
        feature_lbl = conditions.train.features.labels{ff};
        feature = eval(feature_lbl);
        feature_flags(:,ff) = ismember(...
            feature,conditions.train.values.(feature_lbl)(kk,:));
    end
    condition_flags = all(feature_flags,2);
    conditions.train.weights(kk,1) = sum(condition_flags);
end

% normalization
conditions.train.weights = ...
    conditions.train.weights/ max(conditions.train.weights);

% iterate through conditions
for kk = 1 : conditions.test.n
    
    % flag trials for the current condition
    feature_flags = false(n_total_trials,conditions.train.features.n);
    for ff = 1 : conditions.train.features.n
        feature_lbl = conditions.test.features.labels{ff};
        feature = eval(feature_lbl);
        feature_flags(:,ff) = ismember(...
            feature,conditions.test.values.(feature_lbl)(kk,:));
    end
    condition_flags = all(feature_flags,2);
    conditions.test.weights(kk,1) = sum(condition_flags);
end

% normalization
conditions.test.weights = ...
    conditions.test.weights/ max(conditions.test.weights);

%% concatenation settings
n_concats_max = 2^7;
n_concats_train = round(conditions.train.weights * n_concats_max);
n_concats_test = round(conditions.test.weights * n_concats_max);
n_concats_total = sum([n_concats_train; n_concats_test]);

%% time settings
roi = [-500,t_set(end)];
roi_n_bins = range(roi) / psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

%% construct spike rate tensor (time X neurons X concatenations)

% workspace clearance
clear concat_tensor P_tR P_tR_avgs P_tR_avg;

% preallocation
P_tR_avgs = nan(roi_n_bins,roi_n_bins,n_contrasts,n_runs);
map_avgs = nan(roi_n_bins,n_contrasts,n_runs);

concat_stimuli = nan(n_concats_total,n_runs);
concat_contrasts = nan(n_concats_total,n_runs);
concat_choices = nan(n_concats_total,n_runs);
concat_evalset = categorical(nan(n_concats_total,n_runs),[0,1],{'train','test'});

% iterate through runs
for rr = 1 : n_runs
    
    % preallocation
    concat_tensor = nan(roi_n_bins,n_neurons,n_concats_total);
    concat_stimuli = nan(n_concats_total,1);
    concat_contrasts = nan(n_concats_total,1);
    concat_choices = nan(n_concats_total,1);
    concat_evalset = categorical(nan(n_concats_total,1),[0,1],{'train','test'});
    
    % iterate through units
    for nn = 1 : n_neurons
        progressreport(nn,n_neurons,...
            sprintf('sampling concatenations (run %i/%i)',rr,n_runs));
        neuron_flags = data.NeuronNumb == flagged_neurons(nn);
        
        % preallocation
        xval_train_trials_bab = cell(conditions.train.n,conditions.test.n);
        
        % iterate through training conditions
        for kk = 1 : conditions.train.n
            
            % flag trials for the current condition
            feature_flags = false(n_total_trials,conditions.train.features.n);
            for ff = 1 : conditions.train.features.n
                feature_lbl = conditions.train.features.labels{ff};
                feature = eval(feature_lbl);
                feature_flags(:,ff) = ismember(...
                    feature,conditions.train.values.(feature_lbl)(kk,:));
            end
            condition_flags = all(feature_flags,2);
            
            % trial selection
            trial_flags = ...
                valid_flags & ...
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_rates = data.SDF(trial_flags,:);
            
            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi;
            alignment_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + t2(trial_flags);
            chunk_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + roi(2);
            aligned_spkrates = spike_rates';
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [roi_n_bins,n_flagged_trials])';
            
            % detect any test conditions that overlap with the current training condition
            feature_flags = false(conditions.test.n,1);
            for ff = 1 : conditions.test.features.n
                feature_lbl = conditions.test.features.labels{ff};
                feature_flags(:,ff) = any(ismember(...
                    conditions.test.values.(feature_lbl),...
                    conditions.train.values.(feature_lbl)(kk,:)),2);
            end
            xval_condition_flags = all(feature_flags,2);
            
            % handle cross-validation for all detected test conditions
            if any(xval_condition_flags)
                xval_condition_idcs = find(xval_condition_flags)';
                
                % iterate through detected test conditions
                for ii = xval_condition_idcs
                    
                    % detect overlap between training & test trials
                    feature_flags = false(n_flagged_trials,1);
                    for ff = 1 : conditions.test.features.n
                        feature_lbl = conditions.test.features.labels{ff};
                        feature = eval(feature_lbl);
                        feature_flags(:,ff) = any(ismember(...
                            feature(flagged_trials),...
                            conditions.test.values.(feature_lbl)(ii,:)),2);
                    end
                    xval_trial_flags = all(feature_flags,2);
                    
                    % split the conflicting trials into training & test subsets
                    xval_trials = flagged_trials(xval_trial_flags);
                    n_xval_trials = numel(xval_trials);
                    n_train_trials = round(n_xval_trials * 1 / 2);
                    xval_train_idcs = randperm(n_xval_trials,n_train_trials);
                    xval_train_trials_bab{kk,ii} = xval_trials(xval_train_idcs);
                end
                
                % concatenate sub-sampled training sets across test conditions
                train_idcs = find(ismember(...
                    flagged_trials,vertcat(xval_train_trials_bab{kk,:})));
            else
                
                % train using all trials in the remaining conditions
                train_idcs = 1 : n_flagged_trials;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(train_idcs,n_concats_train(kk),true);
            concat_idcs = (1 : n_concats_train(kk)) + ...
                sum(n_concats_train(1:kk-1));
            concat_tensor(:,nn,concat_idcs) = aligned_spkrates(rand_idcs,:)';
            concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
            concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
            concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
            concat_evalset(concat_idcs) = 'train';
        end
        
        % iterate through conditions
        for kk = 1 : conditions.test.n
            
            % flag trials for the current condition
            feature_flags = false(n_total_trials,conditions.test.features.n);
            for ff = 1 : conditions.test.features.n
                feature_lbl = conditions.test.features.labels{ff};
                feature = eval(feature_lbl);
                feature_flags(:,ff) = ismember(...
                    feature,conditions.test.values.(feature_lbl)(kk,:));
            end
            condition_flags = all(feature_flags,2);
            
            % trial selection
            trial_flags = ...
                valid_flags & ...
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_rates = data.SDF(trial_flags,:);
            
            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi;
            alignment_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + t2(trial_flags);
            chunk_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + roi(2);
            aligned_spkrates = spike_rates';
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [roi_n_bins,n_flagged_trials])';
            
            % handle cross-validation
            test_idcs = find(~ismember(...
                flagged_trials,vertcat(xval_train_trials_bab{:,kk})));
            if isempty(test_idcs)
                continue;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(test_idcs,n_concats_test(kk),true);
            concat_idcs = (1 : n_concats_test(kk)) + ...
                sum(n_concats_test(1:kk-1)) + ...
                sum(n_concats_train);
            concat_tensor(:,nn,concat_idcs) = aligned_spkrates(rand_idcs,:)';
            concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
            concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
            concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
            concat_evalset(concat_idcs) = 'test';
        end
    end
    
    %% naive bayes decoder
    nbdopt = struct();
    nbdopt.n_xpoints = 100;
    nbdopt.time = roi_time;
    nbdopt.train.trial_idcs = find(concat_evalset == 'train');
    nbdopt.train.n_trials = numel(nbdopt.train.trial_idcs);
    nbdopt.test.trial_idcs = find(concat_evalset == 'test');
    nbdopt.test.n_trials = numel(nbdopt.test.trial_idcs);
    nbdopt.assumepoissonmdl = false;
    nbdopt.verbose = true;
    
    tic
    [P_tR,~,pthat] = naivebayestimedecoder(concat_tensor,nbdopt);
    toc
    
    %% store current run
    
    % workspace clearance
    clear concat_tensor;
    
    % iterate through contrast conditions
    for ii = 1 : n_contrasts
        contrast_flags = ...
            concat_contrasts(concat_evalset == 'test') == contrast_set(ii);
        P_tR_avgs(:,:,ii,rr) = nanmean(P_tR(:,:,contrast_flags),3);
        map_avgs(:,ii,rr) = nanmean(pthat.mode(:,contrast_flags),2);
    end
end

%% choice of average function
avgfun = @(x,d)nanmean(x,d);
errfun = @(x,d)nanstd(x,0,d);

%% plot superimposed contrast-split posterior averages
fig = figure(...
    figopt,...
    'name',sprintf('superimposed_posterior_averages_%s',contrast_str),...
    'numbertitle','off');
axes(...
    axesopt.default,...
    'xlim',roi2plot,...
    'xtick',unique([roi';roi2plot';0;t_set]),...
    'ylim',roi2plot,...
    'ytick',unique([roi';roi2plot';0;t_set]),...
    'xticklabelrotation',0,...
    'yticklabelrotation',0);
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% convert from tensor to rgb
P_tR_avg = squeeze(avgfun(P_tR_avgs,4));
P = tensor2rgb(permute(P_tR_avg,[2,1,3]),contrast_clrs);
imagesc(roi,roi,P);

% zero lines
plot(xlim,[1,1]*0,':k');

% iterate through slices
for ii = 1 : n_slices
    
    % plot slice handles
    plot([1,1]*slice_times(ii),ylim,':k');
    plot(slice_times(ii),max(ylim),...
        'marker','v',...
        'markersize',7.5,...
        'markerfacecolor',slice_clrs(ii,:),...
        'markeredgecolor','k',...
        'linewidth',1.5);
end

% inset with pseudo colorbar
axes(...
    axesopt.default,...
    'position',[0.075 0.6500 0.2583 0.2717],...
    'yaxislocation','right',...
    'xcolor','none',...
    'xlim',[0,1],...
    'ylim',[0,1],...
    'ytick',0,...
    'colormap',colorlerp(...
    [contrast_clrs(1,:);[1,1,1];contrast_clrs(end,:)],2^8));
ylabel('P(t|R)',...
    'verticalalignment','middle',...
    'rotation',-90);

% colorbar settings
clrbar_width = .05;

% iterate through contrasts
for ii = 1 : n_contrasts
    
    % patch pseudo-colorbar
    xpatch = (1 - clrbar_width * n_contrasts) + ...
        clrbar_width * ((ii - 1) + [0,1,1,0]);
    ypatch = [0,0,1,1];
    patch(xpatch,ypatch,contrast_clrs(ii,:),...
        'edgecolor','none',...
        'linewidth',1.5);
end

% inset with extreme posterior subtractions
axes(...
    axesopt.inset.se,...
    axesopt.default,...
	'xlim',roi2plot,...
    'xtick',unique([roi';roi2plot';0;t_set]),...
    'ylim',roi2plot,...
    'ytick',unique([roi';roi2plot';0;t_set]),...
    'box','on',...
    'colormap',colorlerp(...
    [contrast_clrs(1,:);[1,1,1];contrast_clrs(end,:)],2^8));
title('P(t|R) - P(t|R)');
ylabel('\DeltaP(t|R)',...
    'rotation',-90,...
    'verticalalignment','bottom');

% posterior subtraction
p_contrast_min = avgfun(P_tR(:,:,...
    concat_contrasts(concat_evalset == 'test') == contrast_set(1),:),[3,4]);
p_contrast_max = avgfun(P_tR(:,:,...
    concat_contrasts(concat_evalset == 'test') == contrast_set(end),:),[3,4]);
p_contrast_min = p_contrast_min - min(p_contrast_min,[],[1,2]);
p_contrast_max = p_contrast_max - min(p_contrast_max,[],[1,2]);
p_diff = p_contrast_max - p_contrast_min;
imagesc(roi,roi,p_diff',[-1,1] * n_t / n_tbins * 5);

% zero lines
plot(xlim,[1,1]*0,':k');
plot([1,1]*0,ylim,':k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot superimposed contrast-split MAPs

% figure initialization
fig = figure(...
    figopt,...
    'name','posterior MAPs',...
    'numbertitle','off');

% axes initialization
axes(...
    axesopt.default,...
    'xlim',[-500,t_set(end-2)],...
    'xtick',unique([roi';-500;0;t_set]),...
    'ylim',[-500,t_set(end-2)],...
    'ytick',unique([roi';-500;0;t_set]));
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% iterate through contrast conditions
for ii = 1 : n_contrasts
    map_avg = squeeze(avgfun(map_avgs(:,ii,:),3));
    map_err = squeeze(errfun(map_avgs(:,ii,:),3));
    errorpatch(roi_time,map_avg,map_err,contrast_clrs(ii,:),...
        'facealpha',.25);
    plot(roi_time,map_avg,...
        'color',contrast_clrs(ii,:),...
        'linewidth',1.5);
end

% zero lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*0,':k');

% plot identity line
plot(xlim,ylim,'--k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end