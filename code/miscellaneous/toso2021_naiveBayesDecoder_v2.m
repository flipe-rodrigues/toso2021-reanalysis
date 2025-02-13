%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

% prealocation
P_test = [];

%% iterate through T1 conditions
for tt = 1 : n_t
    
    %% training & test set conditions
    
    % training set conditions
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set(tt),[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
    
    % test set conditions
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),t_set(tt),[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
    
    % print training & test conditions
    fprintf('\nTRAINING CONDITIONS:\n');
    disp(conditions.train.values);
    fprintf('\nTEST CONDITIONS:\n');
    conditions.test.values
    
    %% run settings
    n_runs = 1;
    
    %% concatenation settings
    n_concatspercond = 2^8; % 2^8;
    n_concats = n_concatspercond * (conditions.test.n + conditions.train.n);
    
    %% time settings
    roi = [-0,t_set(end)]; % [-500,t_set(end)] !!!!!!!!!!!
    roi_n_bins = range(roi) / psthbin;
    roi_time = linspace(roi(1),roi(2),roi_n_bins);
    
    %% construct spike rate tensor (time X neurons X concatenations)
    
    % data clearance
    clear concat_tensor P_tR;
    
    % preallocation
    nbd = struct();
    concat_stimuli = nan(n_concats,n_runs);
    concat_contrasts = nan(n_concats,n_runs);
    concat_choices = nan(n_concats,n_runs);
    concat_evalset = categorical(nan(n_concats,n_runs),[0,1],{'train','test'});
    
    % data type selection
    spike_data_field = 'FR';
    
    % iterate through runs
    for rr = 1 : n_runs
        
        % preallocation
        concat_tensor = nan(roi_n_bins,n_neurons,n_concats);
        concat_stimuli = nan(n_concats,1);
        concat_contrasts = nan(n_concats,1);
        concat_choices = nan(n_concats,1);
        concat_evalset = categorical(nan(n_concats,1),[0,1],{'train','test'});
        
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
                spike_counts = data.(spike_data_field)(trial_flags,:);
                spike_rates = ...
                    conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
                
                % S2-offset-aligned spike rates
                alignment = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags) + ...
                    isi;
                alignment_flags = ...
                    valid_time >= alignment + roi(1) & ...
                    valid_time < alignment + t2(trial_flags);
                chunk_flags = ...
                    valid_time >= alignment + roi(1) & ...
                    valid_time < alignment + roi(2);
                aligned_spkrates = spike_rates;
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
                rand_idcs = randsample(train_idcs,n_concatspercond,true);
                concat_idcs = (1 : n_concatspercond) + ...
                    n_concatspercond * (kk - 1);
                concat_tensor(:,nn,concat_idcs) = aligned_spkrates(rand_idcs,:)';
                concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
                concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
                concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
                concat_evalset(concat_idcs) = 'train';
            end
            
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
                spike_counts = data.(spike_data_field)(trial_flags,:);
                spike_rates = ...
                    conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
                
                % S2-offset-aligned spike rates
                alignment = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags) + ...
                    isi;
                alignment_flags = ...
                    valid_time >= alignment + roi(1) & ...
                    valid_time < alignment + t2(trial_flags);
                chunk_flags = ...
                    valid_time >= alignment + roi(1) & ...
                    valid_time < alignment + roi(2);
                aligned_spkrates = spike_rates;
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
                rand_idcs = randsample(test_idcs,n_concatspercond,true);
                concat_idcs = (1 : n_concatspercond) + ...
                    n_concatspercond * (kk + conditions.train.n - 1);
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
        nbdopt.shuffle = false;
        nbdopt.n_shuffles = 10;
        nbdopt.assumepoissonmdl = false;
        nbdopt.verbose = true;
        
        tic
        [P_tR,P_Rt,pthat,neurons] = naivebayestimedecoder(concat_tensor,nbdopt);
        toc
    end
    
    %% choice of average function
    avgfun = @nanmean;
    
    %% plot extreme subtractions of posterior averages
    
    % figure initialization
    fig = figure(...
        figopt,...
        'name','posterior subtractions (extreme)',...
        'numbertitle','off');
    
    % axes initialization
    axes(...
        axesopt.default,...
        'xlim',[0,t_set(end)],...
        'xtick',sort([0;t_set]),...
        'ylim',[0,t_set(end)],...
        'ytick',sort([0;t_set]));
    xlabel('Real time since S_2 onset (ms)');
    ylabel('Decoded time since S_2 onset (ms)');
    
    % posterior subtraction
    stimulus_flags = ismember(...
        concat_stimuli(concat_evalset == 'test'),conditions.test.values.t2);
    contrast_min_flags = ...
        concat_contrasts(concat_evalset == 'test') == contrast_set(1) & ...
        stimulus_flags;
    contrast_max_flags = ...
        concat_contrasts(concat_evalset == 'test') == contrast_set(end) & ...
        stimulus_flags;
    p_contrast_min = avgfun(P_tR(:,:,contrast_min_flags),3);
    p_contrast_min = p_contrast_min ./ nansum(p_contrast_min,2);
    p_contrast_max = avgfun(P_tR(:,:,contrast_max_flags),3);
    p_contrast_max = p_contrast_max ./ nansum(p_contrast_max,2);
    p_diff = p_contrast_max - p_contrast_min;
    imagesc([0,t_set(end)],[0,t_set(end)],p_diff',...
        [-1,1] * n_t / n_tbins * 1);
    plot(xlim,ylim,'--w');
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
    
    %% plot contrast-split posterior averages
    figure(...
        'name','contrast-split posterior averages',...
        'numbertitle','off',...
        'windowstyle','docked');
    sps = gobjects(n_contrasts,1);
    for ii = 1 : n_contrasts
        sps(ii) = subplot(1,n_contrasts,ii);
    end
    set(sps,...
        axesopt.default,...
        'xlim',[0,t_set(end)],...
        'xtick',sort([0;t_set]),...
        'ylim',[0,t_set(end)],...
        'ytick',sort([0;t_set]));
    xlabel(sps(contrast_mode_idx),'Real time since S_2 onset (ms)');
    ylabel(sps(1),'Decoded time since S_2 onset (ms)');
    linkaxes(sps);
    
    % iterate through contrast conditions
    for ii = 1 : n_contrasts
        contrast_flags = concat_contrasts(concat_evalset == 'test') == contrast_set(ii);
        concat_flags = ...
            contrast_flags;
        if sum(concat_flags) == 0
            continue;
        end
        title(sps(ii),sprintf('%s = %.0f mm/s',...
            contrast_lbl,contrast_set(ii)),...
            'fontsize',10);
        p_cond = avgfun(P_tR(:,:,concat_flags),3);
        p_cond = p_cond ./ nansum(p_cond,2);
        p_cond(isnan(p_cond)) = max(p_cond(:));
        imagesc(sps(ii),[0,t_set(end)],[0,t_set(end)],p_cond');
        plot(sps(ii),xlim,ylim,'--w');
    end
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
    
    %% plot contrast-split point estimate averages
    pthat_avgfuns = {...
        ...@(x)nanmean(x,2),...
        @(x)nanmedian(x,2),...
        };
    pthat_errfuns = {...
        ...@(x)nanmean(x,2) + [-1,1] .* nanstd(x,0,2) / sqrt(size(x,2));
        @(x)quantile(x,[.25,.75],2),...
        };

    % fade settings
    fadeifnoisy = true;
    alphabounds_sem = [.05,.3];
    alphabounds_mu = [.15,.85];
    
    % point estimate selection
    type = 'mode';
    
    % figure initialization
    fig = figure(...
        figopt,...
        'name',sprintf('point estimates (%s)',type),...
        'numbertitle','off');
    
    % axes initialization
    axes(...
        axesopt.default,...
        'layer','bottom',...
        'xlim',[0,t_set(end)],...
        'xtick',sort([0;t_set]),...
        'ylim',[0,t_set(end)],...
        'ytick',sort([0;t_set]));
    xlabel('Real time since S_2 onset (ms)');
    ylabel('Decoded time since S_2 onset (ms)');
    
    % plot unity line
    plot(xlim,ylim,'--k');
    
    % iterate through contrast conditions
    for ii = 1 : n_contrasts
        contrast_flags = concat_contrasts(concat_evalset == 'test') == contrast_set(ii);
        stimulus_flags = ...
            ismember(concat_stimuli(concat_evalset == 'test'),conditions.test.values.t2);
        concat_flags = ...
            contrast_flags & ...
            stimulus_flags;
        if sum(concat_flags) == 0
            continue;
        end
        
        % compute point estimate average & error
        pthat_avg = pthat_avgfuns{1}(pthat.(type)(:,concat_flags));
        pthat_err = pthat_errfuns{1}(pthat.(type)(:,concat_flags));
        
        % compute surviving trial count
        if fadeifnoisy
            trials_throughtime = sum(~isnan(pthat.(type)),2)';
        end
        
        % patch error
        xpatch = [roi_time,fliplr(roi_time)];
        ypatch = [pthat_err(:,1);flipud(pthat_err(:,2))];
        if fadeifnoisy
            patch_alpha = [trials_throughtime,fliplr(trials_throughtime)];
            patch_alpha = (patch_alpha - min(trials_throughtime)) / ...
                range(trials_throughtime);
            patch_alpha = patch_alpha * alphabounds_sem(2) + alphabounds_sem(1);
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                patch(...
                    xpatch(alpha_flags),...
                    ypatch(alpha_flags),0,...
                    'facealpha',alpha_levels(aa),...
                    'edgecolor','none',...
                    'facecolor',contrast_clrs(ii,:));
            end
        else
            patch(xpatch,ypatch,contrast_clrs(ii,:),...
                'facecolor',contrast_clrs(ii,:),...
                'edgecolor','none',...
                'facealpha',.25);
        end
        
        % patch average
        if fadeifnoisy
            patch_alpha = trials_throughtime;
            patch_alpha = (patch_alpha - min(trials_throughtime)) / ...
                range(trials_throughtime);
            patch_alpha = patch_alpha * alphabounds_mu(2) + alphabounds_mu(1);
            alpha_levels = unique(patch_alpha,'stable');
            n_alpha_levels = numel(alpha_levels);
            for aa = 1 : n_alpha_levels
                alpha_flags = patch_alpha == alpha_levels(aa);
                patch(...
                    [roi_time(alpha_flags),nan],...
                    [pthat_avg(alpha_flags)',nan],0,...
                    'edgealpha',alpha_levels(aa),...
                    'edgecolor',contrast_clrs(ii,:),...
                    'facecolor','none',...
                    'linewidth',1.5);
            end
        else
            plot(roi_time,pthat_avg,...
                'color',contrast_clrs(ii,:),...
                'linewidth',1.5);
        end
        
        % plot stimulus onset
        onset_idx = 1;
        plot(roi_time(onset_idx),pthat_avg(onset_idx),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',7.5,...
            'markerfacecolor','w',...
            'markeredgecolor',contrast_clrs(ii,:));
        
        % iterate through stimuli
%         for jj = 1 : n_t
%             
%             % plot stimulus offset
%             offset_idx = find(roi_time >= t_set(jj),1) - 1;
%             scatter(roi_time(offset_idx),pthat_avg(offset_idx),60,...
%                 'linewidth',1.5,...
%                 'marker','o',...
%                 'markerfacealpha',alpha_levels(jj).^fadeifnoisy,...
%                 'markerfacecolor',contrast_clrs(ii,:),...
%                 'markeredgecolor','none');
%         end
    end
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end