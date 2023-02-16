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
    
    % save current run
    %     nbd.P_tR(:,:,:,rr) = P_tR;
    %     nbd.P_Rt(:,:,:,rr) = P_Rt;
    %     nbd.pthat.mode(:,:,rr) = pthat.mode;
    %     nbd.pthat.median(:,:,rr) = pthat.median;
    %     nbd.pthat.mean(:,:,rr) = pthat.mean;
    %     nbd.neurons(:,rr) = neurons;
end

%% plot likelihoods
figure;
set(gca,...
    axesopt.default,...,...
    'xlim',[0,t_set(end)],...
    'xtick',sort([0;t_set]));
xlabel('Time since T_2 (ms)');
ylabel('Firing rate (Hz)');

% iterate through units
for nn = 1 : n_neurons
    cla;
    r_bounds = neurons(nn).x_bounds;
    r_bw = neurons(nn).x_bw;
    if range(r_bounds) == 0
        continue;
    end
    ylim(r_bounds);
    title(sprintf('neuron: %i, bw: %.2f',nn,r_bw));
    p_Rt = squeeze(P_Rt(:,nn,:));
    nan_flags = isnan(p_Rt);
    if sum(abs(diff(any(nan_flags,2)))) > 1
        fprintf('\tcheck neuron %i!\n',nn);
    end
    p_Rt(nan_flags) = max(p_Rt(:));
    imagesc(xlim,r_bounds,p_Rt');
    for ii = 1 : n_t
        plot([1,1]*t_set(ii),ylim,'--w');
    end
    drawnow;
    pause(.1);
end

%% plot single-trial posteriors
figure;
set(gca,...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',sort([0;t_set]));
xlabel('Real time (ms)');
ylabel('Decoded time (ms)');

% iterate through trials
for kk = randperm(nbdopt.test.n_trials,min(nbdopt.test.n_trials,100))
    cla;
    title(sprintf('trial: %i, T_2: %i, T_2: %i',kk,...
        concat_stimuli(nbdopt.test.trial_idcs(kk)),...
        concat_contrasts(nbdopt.test.trial_idcs(kk))));
    p_tR = squeeze(P_tR(:,:,kk));
    p_tR(isnan(p_tR)) = max(p_tR(:));
    imagesc(xlim,ylim,p_tR');
    plot(roi_time,pthat.mode(:,kk),...
        'color','w',...
        'linestyle','-',...
        'linewidth',1);
    plot(roi_time,pthat.median(:,kk),...
        'color','c',...
        'linestyle','-',...
        'linewidth',1);
    plot(roi_time,pthat.mean(:,kk),...
        'color','r',...
        'linestyle','-',...
        'linewidth',1);
    for ii = 1 : n_t
        plot([1,1]*t_set(ii),ylim,'--w');
    end
    pause(.1);
    drawnow;
end

%% choice of average function
avgfun = @nanmean;

%% plot contrast- & stimulus-split posterior averages
figure(...
    'name','stimulus-split posterior averages',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(n_contrasts,n_stimuli);
for ii = 1 : n_contrasts
    for jj = 1 : n_stimuli
        sp_idx = jj + (ii - 1) * n_stimuli;
        sps(ii,jj) = subplot(n_contrasts,n_stimuli,sp_idx);
        %         xlabel(sps(ii,jj),'Real time since S_2 onset (ms)');
        %         ylabel(sps(ii,jj),'Decoded time since S_2 onset (ms)');
    end
end
set(sps,...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',[],...sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',[]); % sort([0;t_set]));
linkaxes(sps);

% iterate through contrast conditions
for ii = 1 : n_contrasts
    contrast_flags = concat_contrasts(concat_evalset == 'test') == contrast_set(ii);
    
    % iterate through stimuli
    for jj = 1 : n_stimuli
        stimulus_flags = concat_stimuli(concat_evalset == 'test') == stim_set(jj);
        concat_flags = ...
            contrast_flags & ...
            stimulus_flags;
        if sum(concat_flags) == 0
            continue;
        end
        title(sps(ii,jj),sprintf('T_2 = %.0f ms, %s = %.0f mm/s',...
            t_set(jj),contrast_lbl,contrast_set(ii)),...
            'fontsize',10);
        p_cond = avgfun(P_tR(:,:,concat_flags),3);
        p_cond = p_cond ./ nansum(p_cond,2);
        p_cond(isnan(p_cond)) = max(p_cond(:));
        imagesc(sps(ii,jj),[0,t_set(end)],[0,t_set(end)],p_cond');
        plot(sps(ii,jj),xlim,ylim,'--w');
    end
end

%% plot intermediate subtractions of stimulus-split posterior averages
figure(...
    'name','stimulus-split posterior subtractions (intermediate)',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(n_contrasts,n_stimuli);
for ii = 1 : n_contrasts
    for jj = 1 : n_stimuli
        sp_idx = jj + (ii - 1) * n_stimuli;
        sps(ii,jj) = subplot(n_contrasts,n_stimuli,sp_idx);
        %         xlabel(sps(ii,jj),'Real time since S_2 onset (ms)');
        %         ylabel(sps(ii,jj),'Decoded time since S_2 onset (ms)');
    end
end
set(sps,...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',[],...sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',[]); % sort([0;t_set]));
linkaxes(sps);

% flag reference intensity concatenations
contrast_ref_flags = concat_contrasts(concat_evalset == 'test') == contrast_set(contrast_mode_idx);

% iterate through contrast conditions
for ii = 1 : n_contrasts
    contrast_flags = concat_contrasts(concat_evalset == 'test') == contrast_set(ii);
    
    % iterate through stimuli
    for jj = 1 : n_stimuli
        stimuli_flags = concat_stimuli(concat_evalset == 'test') == stim_set(jj);
        ref_flags = ...
            contrast_ref_flags & ...
            stimuli_flags;
        concat_flags = ...
            contrast_flags & ...
            stimuli_flags;
        if sum(ref_flags) == 0 || sum(concat_flags) == 0
            continue;
        end
        title(sps(ii,jj),sprintf('T_2: %.0f ms, %s: %.0f - %.0f mm/s',...
            t_set(jj),contrast_lbl,...
            contrast_set(ii),contrast_set(contrast_mode_idx)),...
            'fontsize',10);
        p_ref = avgfun(P_tR(:,:,ref_flags),3);
        p_ref = p_ref ./ nansum(p_ref,2);
        p_cond = avgfun(P_tR(:,:,concat_flags),3);
        p_cond = p_cond ./ nansum(p_cond,2);
        p_diff = p_cond - p_ref;
        imagesc(sps(ii,jj),[0,t_set(end)],[0,t_set(end)],p_diff',...
            [-1,1] * n_t / n_tbins * 1);
        plot(sps(ii,jj),xlim,ylim,'--w');
    end
end

%% plot extreme subtractions of stimulus-split posterior averages
figure(...
    'name','stimulus-split posterior subtractions (extreme)',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(n_stimuli,1);
for jj = 1 : n_stimuli
    sps(jj) = subplot(1,n_stimuli,jj);
    xlabel(sps(jj),'Real time since S_2 onset (ms)');
end
ylabel(sps(1),'Decoded time since S_2 onset (ms)');
set(sps,...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',sort([0;t_set]));
linkaxes(sps);

% iterate through stimuli
for jj = 1 : n_stimuli
    stimulus_flags = concat_stimuli(concat_evalset == 'test') == stim_set(jj);
    contrast_min_flags = ...
        concat_contrasts(concat_evalset == 'test') == contrast_set(1) & ...
        stimulus_flags;
    contrast_max_flags = ...
        concat_contrasts(concat_evalset == 'test') == contrast_set(end) & ...
        stimulus_flags;
    title(sps(jj),sprintf('%s: %.0f - %.0f mm/s',...
        contrast_lbl,contrast_set(end),contrast_set(1)),...
        'fontsize',10);
    p_contrast_min = avgfun(P_tR(:,:,contrast_min_flags),3);
    p_contrast_min = p_contrast_min ./ nansum(p_contrast_min,2);
    p_contrast_max = avgfun(P_tR(:,:,contrast_max_flags),3);
    p_contrast_max = p_contrast_max ./ nansum(p_contrast_max,2);
    p_diff = p_contrast_max - p_contrast_min;
    imagesc(sps(jj),[0,t_set(end)],[0,t_set(end)],p_diff',...
        [-1,1] * n_t / n_tbins * 1);
    plot(sps(jj),xlim,ylim,'--w');
end

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
    'xtick',[],...sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',[]); % sort([0;t_set]));
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

%% plot superimposed contrast-split posterior averages
figure(...
    'name','posterior averages',...
    'numbertitle','off');
axes(...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',sort([0;t_set]),...
    'colormap',colorlerp([[1,1,1];contrast_clrs(ii,:)],2^8));
xlabel('Real time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');
clims = [0,.05]; % quantile(avgfun(P_tR,3),[0,1],[1,2])';

% iterate through contrast conditions
for ii = 1 : n_contrasts
    contrast_flags = concat_contrasts(concat_evalset == 'test') == contrast_set(ii);
    concat_flags = ...
        contrast_flags;
    if sum(concat_flags) == 0
        continue;
    end
    p_cond = avgfun(P_tR(:,:,concat_flags),3);
    p_cond = p_cond ./ nansum(p_cond,2);
    p_cond(p_cond < clims(1)) = clims(1);
    p_cond(p_cond > clims(2)) = clims(2);
    p_cond(isnan(p_cond)) = max(p_cond(:));
    p_patch = mat2patch(p_cond,[0,t_set(end)],[0,t_set(end)]);
    patch(p_patch,...
        'facevertexalphadata',p_patch.facevertexcdata,...
        'edgecolor','none',...
        'facecolor',contrast_clrs(ii,:),...
        'facealpha','flat');
end

% plot identity line
plot(xlim,ylim,'--k');

%% plot contrast- & stimulus-split point estimate averages
pthat_avgfuns = {...
    ...@(x)nanmean(x,2),...
    @(x)nanmedian(x,2),...
    };
pthat_errfuns = {...
    ...@(x)nanmean(x,2) + [-1,1] .* nanstd(x,0,2) / sqrt(size(x,2));
    @(x)quantile(x,[.25,.75],2),...
    };
n_funs = numel(pthat_avgfuns);

% iterate through point estimate types
pthat_types = fieldnames(pthat);
n_pthats = numel(pthat_types);
for tt = 1 : n_pthats
    type = pthat_types{tt};
    
    % figure initialization
    figure(...
        'name',sprintf('stimulus-split point estimates (%s)',type),...
        'numbertitle','off',...
        'windowstyle','docked');
    n_rows = n_funs * 2;
    n_cols = ceil(n_stimuli / 2);
    n_sps = n_rows * n_cols;
    sps = gobjects(n_sps,1);
    for ii = 1 : n_rows
        for jj = 1 : n_cols
            sp_idx = jj + (ii - 1) * n_cols;
            sps(sp_idx) = subplot(n_rows,n_cols,sp_idx);
            xlabel(sps(sp_idx),'Real time since S_2 onset (ms)');
        end
        ylabel(sps(sp_idx),'Decoded time since S_2 onset (ms)');
    end
    set(sps,...
        axesopt.default,...
        'xlim',[0,t_set(end)],...
        'xtick',sort([0;t_set]),...
        'ylim',[0,t_set(end)],...
        'ytick',sort([0;t_set]));
    linkaxes(sps);
    
    % iterate through average & error functions
    for ff = 1 : n_funs
        
        % iterate through contrast conditions
        for ii = 1 : n_contrasts
            contrast_flags = concat_contrasts(concat_evalset == 'test') == contrast_set(ii);
            
            % iterate through stimuli
            for jj = 1 : n_stimuli
                sp_idx = jj + (ff - 1) * n_cols;
                stimulus_flags = concat_stimuli(concat_evalset == 'test') == stim_set(jj);
                concat_flags = ...
                    contrast_flags & ...
                    stimulus_flags;
                if sum(concat_flags) == 0
                    continue;
                end
                
                title(sps(sp_idx),sprintf('T_2 = %.2f ms',t_set(jj)));
                
                pthat_avg = pthat_avgfuns{ff}(pthat.(type)(:,concat_flags));
                pthat_err = pthat_errfuns{ff}(pthat.(type)(:,concat_flags));
                
                t2bin_idcs = 1 : t_set(jj) / psthbin;
                
                xpatch = [roi_time(t2bin_idcs),fliplr(roi_time(t2bin_idcs))];
                ypatch = [pthat_err(t2bin_idcs,1);flipud(pthat_err(t2bin_idcs,2))];
                patch(sps(sp_idx),...
                    xpatch,ypatch,contrast_clrs(ii,:),...
                    'facecolor',contrast_clrs(ii,:),...
                    'edgecolor','none',...
                    'facealpha',.25);
                plot(sps(sp_idx),...
                    roi_time(t2bin_idcs),pthat_avg(t2bin_idcs),...
                    'color',contrast_clrs(ii,:),...
                    'linewidth',1.5);
                
                onset_idx = 1;
                plot(sps(sp_idx),...
                    roi_time(onset_idx),pthat_avg(onset_idx),...
                    'color',contrast_clrs(ii,:),...
                    'linewidth',1.5,...
                    'marker','o',...
                    'markersize',8.5,...
                    'markerfacecolor','w',...
                    'markeredgecolor',contrast_clrs(ii,:));
                offset_idx = find(roi_time >= t_set(jj),1) - 1;
                plot(sps(sp_idx),...
                    roi_time(offset_idx),pthat_avg(offset_idx),...
                    'color',contrast_clrs(ii,:),...
                    'linewidth',1.5,...
                    'marker','o',...
                    'markersize',8.5,...
                    'markerfacecolor',contrast_clrs(ii,:),...
                    'markeredgecolor','k');
                
                % identity line
                plot(sps(sp_idx),xlim,ylim,'--k');
            end
        end
    end
end

%% plot contrast-split point estimate averages

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
    pthat_avg = pthat_avgfuns{ff}(pthat.(type)(:,concat_flags));
    pthat_err = pthat_errfuns{ff}(pthat.(type)(:,concat_flags));
    
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
    for jj = 1 : n_t

        % plot stimulus offset
        offset_idx = find(roi_time >= t_set(jj),1) - 1;
        scatter(roi_time(offset_idx),pthat_avg(offset_idx),60,...
            'linewidth',1.5,...
            'marker','o',...
            'markerfacealpha',alpha_levels(jj).^fadeifnoisy,...
            'markerfacecolor',contrast_clrs(ii,:),...
            'markeredgecolor','none');
    end
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end