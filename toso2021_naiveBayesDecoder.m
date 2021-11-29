%% initialization
if ~exist('tosoData','var')
    toso2021_wrapper;
end

%% stimulus selection (for training & test sets)

% stimuli
stim2train_idcs = [1 : 2, 6 : n_stimuli];
stim2train_n = numel(stim2train_idcs);
stim2test_idcs = [3,4,5];
stim2test_n = numel(stim2test_idcs);

% contrasts
contrast2train_idcs = contrast_mode_idx;
contrast2train_n = numel(contrast2train_idcs);
contrast2test_idcs = 1 : n_contrasts;
contrast2test_n = numel(contrast2test_idcs);

%% condition properties (for training & test sets)
conditions = struct();

% training set conditions
conditions.train.n = stim2train_n * contrast2train_n;
conditions.train.idcs = ...
    num2cell(repmat(stim2train_idcs',contrast2train_n,1));
conditions.train.stimuli = cellfun(...
    @(x) stim_set(x),conditions.train.idcs,...
    'uniformoutput',false);
conditions.train.contrasts = ...
    num2cell(repmat(contrast_set(contrast2train_idcs'),stim2train_n,1));

% test set conditions
conditions.test.n = stim2test_n * contrast2test_n;
conditions.test.idcs = ...
    num2cell(sort(repmat(stim2test_idcs',contrast2test_n,1)));
conditions.test.stimuli = cellfun(...
    @(x) stim_set(x),conditions.test.idcs,...
    'uniformoutput',false);
conditions.test.contrasts = ...
    num2cell(repmat(contrast_set(contrast2test_idcs'),stim2test_n,1));

%% concatenation settings
n_concatspercond = 2^8;
n_concats = n_concatspercond * (conditions.test.n + conditions.train.n);

%% construct spike rate tensor (time X neurons X concatenations)
roi = [0,t_set(end)];
roi_n_bins = range(roi) * psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

% preallocation
concat_tensor = nan(roi_n_bins,n_neurons,n_concats);
concat_contrasts = nan(n_concats,1);
concat_stimuli = nan(n_concats,1);
concat_choices = nan(n_concats,1);
concat_evalset = nan(n_concats,1);

% stimulus 1 flags
i1_flags = i1 == i_set(i1_mode_idx);
t1_flags = t1 > t_set(2) & t1 < t_set(end-1);

% iterate through units
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing concatenations');
    neuron_flags = tosoData.NeuronNumb == flagged_neurons(nn);
    
    % preallocation
    nn_train_trials = cell(conditions.train.n,1);
    
    % iterate through training conditions
    for kk = 1 : conditions.train.n
        
        % flag trials for the current condition
        stimulus_flags = ismember(stimuli,conditions.train.stimuli{kk});
        contrast_flags = ismember(contrasts,conditions.train.contrasts{kk});
        s2_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            stimulus_flags & ...
            ...pseudo_correct_flags & ...
            i1_flags & ...
            t1_flags & ...
            contrast_flags;
        flagged_trials = find(s2_spike_flags);
        if sum(s2_spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        s2_spike_counts = tosoData.FR(s2_spike_flags,:);
        s2_spike_rates = ...
            conv2(1,t_kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(s2_spike_counts,1);
        
        % T2-offset-aligned spike rates
        s2_alignment_offset = ...
            pre_init_padding + ...
            pre_t1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            inter_t1t2_delay;
        s2_alignment_flags = ...
            t_valid_time >= s2_alignment_offset + roi(1) & ...
            t_valid_time < s2_alignment_offset + t2(s2_spike_flags);
        s2_chunk_flags = ...
            t_valid_time >= s2_alignment_offset + roi(1) & ...
            t_valid_time < s2_alignment_offset + roi(2);
        s2_spkrates = s2_spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = ...
            reshape(s2_spkrates(s2_chunk_flags'),[roi_n_bins,n_trials])';
        
        % handle cross-validation
        if ismember(conditions.train.stimuli{kk},[conditions.test.stimuli{:}]) && ...
                ismember(conditions.train.contrasts{kk},[conditions.test.contrasts{:}])
            n_train_trials = round(n_trials * 2 / 3);
            train_idcs = randperm(n_trials,n_train_trials);
        else
            n_train_trials = n_trials;
            train_idcs = 1 : n_train_trials;
        end
        nn_train_trials{kk} = flagged_trials(train_idcs);
        
        % store tensor & concatenation data
        rand_idcs = randsample(train_idcs,n_concatspercond,true);
        concat_idcs = (1 : n_concatspercond) + ...
            n_concatspercond * (kk - 1);
        concat_tensor(:,nn,concat_idcs) = ...
            s2_spkrates(rand_idcs,:)';
        concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
        concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
        concat_choices(concat_idcs) = choices(flagged_trials(rand_idcs));
        concat_evalset(concat_idcs) = 0;
    end
    
    % vectorize cell with training trial indices
    nn_train_trials = vertcat(nn_train_trials{:});
    
    % iterate through conditions
    for kk = 1 : conditions.test.n
        
        % flag trials for the current condition
        stimulus_flags = ismember(stimuli,conditions.test.stimuli{kk});
        contrast_flags = ismember(contrasts,conditions.test.contrasts{kk});
        s2_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            stimulus_flags & ...
            contrast_flags;
        flagged_trials = find(s2_spike_flags);
        if sum(s2_spike_flags) == 0
            disp([conditions.test.stimuli{kk},...
                conditions.test.contrasts{kk}])
            continue;
        end
        
        % fetch spike counts & compute spike rates
        s2_spike_counts = tosoData.FR(s2_spike_flags,:);
        s2_spike_rates = ...
            conv2(1,t_kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(s2_spike_counts,1);
        
        % T2-aligned spike rates
        s2_alignment_offset = ...
            pre_init_padding + ...
            pre_t1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            inter_t1t2_delay;
        s2_alignment_flags = ...
            t_valid_time >= s2_alignment_offset + roi(1) & ...
            t_valid_time < s2_alignment_offset + t2(s2_spike_flags);
        s2_chunk_flags = ...
            t_valid_time >= s2_alignment_offset + roi(1) & ...
            t_valid_time < s2_alignment_offset + roi(2);
        s2_spkrates = s2_spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = ...
            reshape(s2_spkrates(s2_chunk_flags'),[roi_n_bins,n_trials])';
        
        % handle cross-validation
        test_flags = ~ismember(flagged_trials,nn_train_trials);
        test_idcs = find(test_flags);
        
        % store tensor & concatenation data
        rand_idcs = randsample(test_idcs,n_concatspercond,true);
        concat_idcs = (1 : n_concatspercond) + ...
            n_concatspercond * (kk + conditions.train.n - 1);
        concat_tensor(:,nn,concat_idcs) = ...
            s2_spkrates(rand_idcs,:)';
        concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
        concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
        concat_choices(concat_idcs) = choices(flagged_trials(rand_idcs));
        concat_evalset(concat_idcs) = 1;
    end
end

% convert train / test flags to categorical
concat_evalset = categorical(concat_evalset,[0,1],{'train','test'});

%% naive bayes decoder
nbdopt = struct();
nbdopt.n_xpoints = 100;
nbdopt.condition = concat_contrasts;
nbdopt.stimulus = concat_stimuli;
nbdopt.time = roi_time;
nbdopt.train.trial_idcs = find(concat_evalset == 'train');
nbdopt.train.n_trials = numel(nbdopt.train.trial_idcs);
nbdopt.test.trial_idcs = find(concat_evalset == 'test');
nbdopt.test.n_trials = numel(nbdopt.test.trial_idcs);
nbdopt.assumepoissonmdl = false;

tic
[P_tR,P_Rt,pthat,neurons] = ...
    naivebayesdecoder_simpler(concat_tensor,nbdopt);
n_neurons = numel(neurons);
toc

%% plot likelihoods
figure;
set(gca,...
    axesopt.default,...,...
    'xlim',[0,t_set(end)],...
    'xtick',sort([0;t_set]));
xlabel('Time since T_2 (ms)');
ylabel('Firing rate (Hz)');

% iterate through units
for nn = 87%1 : n_neurons
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
    plot(nbdopt.time,pthat.mode(:,kk),...
        'color','w',...
        'linestyle','-',...
        'linewidth',1);
    plot(nbdopt.time,pthat.median(:,kk),...
        'color','c',...
        'linestyle','-',...
        'linewidth',1);
    plot(nbdopt.time,pthat.mean(:,kk),...
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

%% plot condition-split posterior averages
figure(...
    'name','posterior averages',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(n_i,stim2test_n);
for ii = 1 : n_i
    for jj = 1 : stim2test_n
        sp_idx = jj + (ii - 1) * stim2test_n;
        sps(ii,jj) = subplot(n_i,stim2test_n,sp_idx);
        %         xlabel(sps(ii,jj),'Real time since T_2 onset (ms)');
        %         ylabel(sps(ii,jj),'Decoded time since T_2 onset (ms)');
    end
end
set(sps,...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',[],...sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',[]); % sort([0;t_set]));
linkaxes(sps);

% iterate through intensities
for ii = 1 : n_i
    i2_flags = concat_contrasts(nbdopt.test.trial_idcs) == i_set(ii);
    
    % iterate through stimuli
    for jj = 1 : stim2test_n
        t2_flags = concat_stimuli(nbdopt.test.trial_idcs) == t_set(stim2test_idcs(jj));
        concat_flags = ...
            t2_flags & ...
            i2_flags;
        title(sps(ii,jj),sprintf('T_2 = %.2f ms, I_2 = %.2f mm/s',...
            t_set(stim2test_idcs(jj)),i_set(ii)));
        p_cond = avgfun(P_tR(:,:,concat_flags),3);
        p_cond = p_cond ./ nansum(p_cond,2);
        p_cond(isnan(p_cond)) = max(p_cond(:));
        imagesc(sps(ii,jj),[0,t_set(end)],[0,t_set(end)],p_cond');
        plot(sps(ii,jj),xlim,ylim,'--w');
    end
end

%% plot control-subtracted posterior averages
figure(...
    'name','posterior subtractions',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(n_i,stim2test_n);
for ii = 1 : n_i
    for jj = 1 : stim2test_n
        sp_idx = jj + (ii - 1) * stim2test_n;
        sps(ii,jj) = subplot(n_i,stim2test_n,sp_idx);
        %         xlabel(sps(ii,jj),'Real time since T_2 onset (ms)');
        %         ylabel(sps(ii,jj),'Decoded time since T_2 onset (ms)');
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
i2_ref_flags = concat_contrasts(nbdopt.test.trial_idcs) == i_set(contrast2train_idcs);

% iterate through intensities
for ii = 1 : n_i
    i2_flags = concat_contrasts(nbdopt.test.trial_idcs) == i_set(ii);
    
    % iterate through stimuli
    for jj = 1 : stim2test_n
        t2_flags = concat_stimuli(nbdopt.test.trial_idcs) == t_set(stim2test_idcs(jj));
        ref_flags = ...
            t2_flags & ...
            i2_ref_flags;
        concat_flags = ...
            t2_flags & ...
            i2_flags;
        
        title(sps(ii,jj),sprintf('T_2: %.2f ms, I_2: %.2f - %.2f mm/s',...
            t_set(stim2test_idcs(jj)),i_set(ii),i_set(contrast2train_idcs)));
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

%% plot subtractions of extreme posterior averages
figure(...
    'name','posterior subtractions',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(stim2test_n,1);
for jj = 1 : stim2test_n
    sps(jj) = subplot(1,stim2test_n,jj);
    xlabel(sps(jj),'Real time since T_2 onset (ms)');
end
ylabel(sps(1),'Decoded time since T_2 onset (ms)');
set(sps,...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',sort([0;t_set]));
linkaxes(sps);

% iterate through stimuli
for jj = 1 : stim2test_n
    t2_flags = concat_stimuli(nbdopt.test.trial_idcs) == t_set(stim2test_idcs(jj));
    i2_min_flags = ...
        t2_flags & ...
        concat_contrasts(nbdopt.test.trial_idcs) == i_set(1);
    i2_max_flags = ...
        t2_flags & ...
        concat_contrasts(nbdopt.test.trial_idcs) == i_set(end);
    title(sps(jj),sprintf('I_2: %.2f - %.2f mm/s',i_set(end),i_set(1)));
    p_i2_min = avgfun(P_tR(:,:,i2_min_flags),3);
    p_i2_min = p_i2_min ./ nansum(p_i2_min,2);
    p_i2_max = avgfun(P_tR(:,:,i2_max_flags),3);
    p_i2_max = p_i2_max ./ nansum(p_i2_max,2);
    p_diff = p_i2_max - p_i2_min;
    imagesc(sps(jj),[0,t_set(end)],[0,t_set(end)],p_diff',...
        [-1,1] * n_t / n_tbins * 1);
    plot(sps(jj),xlim,ylim,'--w');
end

%% plot I2-split point estimate averages
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
        'name',sprintf('point estimates (%s)',type),...
        'numbertitle','off',...
        'windowstyle','docked');
    n_rows = n_funs * 2;
    n_cols = ceil(stim2test_n / 2);
    n_sps = n_rows * n_cols;
    sps = gobjects(n_sps,1);
    for ii = 1 : n_rows
        for jj = 1 : n_cols
            sp_idx = jj + (ii - 1) * n_cols;
            sps(sp_idx) = subplot(n_rows,n_cols,sp_idx);
            xlabel(sps(sp_idx),'Real time since T_2 onset (ms)');
        end
        ylabel(sps(sp_idx),'Decoded time since T_2 onset (ms)');
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
        
        % iterate through intensities
        for ii = 1 : n_i
            i2_flags = concat_contrasts(nbdopt.test.trial_idcs) == i_set(ii);
            
            % iterate through stimuli
            for jj = 1 : stim2test_n
                sp_idx = jj + (ff - 1) * n_cols;
                t2_flags = concat_stimuli(nbdopt.test.trial_idcs) == t_set(stim2test_idcs(jj));
                concat_flags = ...
                    t2_flags & ...
                    i2_flags;
                title(sps(sp_idx),sprintf('T_2 = %.2f ms',t_set(stim2test_idcs(jj))));
                
                pthat_avg = pthat_avgfuns{ff}(pthat.(type)(:,concat_flags));
                pthat_err = pthat_errfuns{ff}(pthat.(type)(:,concat_flags));
                
                t2bin_idcs = (1 : t_set(stim2test_idcs(jj))) * psthbin;
                
                xpatch = [nbdopt.time(t2bin_idcs),fliplr(nbdopt.time(t2bin_idcs))];
                ypatch = [pthat_err(t2bin_idcs,1);flipud(pthat_err(t2bin_idcs,2))];
                patch(sps(sp_idx),...
                    xpatch,ypatch,contrast_clrs(ii,:),...
                    'facecolor',contrast_clrs(ii,:),...
                    'edgecolor','none',...
                    'facealpha',.25+(i_set(ii) == i_set(end))*.15);
                plot(sps(sp_idx),...
                    nbdopt.time(t2bin_idcs),pthat_avg(t2bin_idcs),...
                    'color',contrast_clrs(ii,:),...
                    'linewidth',1.5);
                
                onset_idx = 1;
                plot(sps(sp_idx),...
                    nbdopt.time(onset_idx),pthat_avg(onset_idx),...
                    'color',contrast_clrs(ii,:),...
                    'linewidth',1.5,...
                    'marker','o',...
                    'markersize',8.5,...
                    'markerfacecolor','w',...
                    'markeredgecolor',contrast_clrs(ii,:));
                offset_idx = find(nbdopt.time >= t_set(stim2test_idcs(jj)),1) - 1;
                plot(sps(sp_idx),...
                    nbdopt.time(offset_idx),pthat_avg(offset_idx),...
                    'color',contrast_clrs(ii,:),...
                    'linewidth',1.5,...
                    'marker','o',...
                    'markersize',8.5,...
                    'markerfacecolor',contrast_clrs(ii,:),...
                    'markeredgecolor','k');
                
                plot(sps(sp_idx),xlim,ylim,'--k');
            end
        end
    end
end

%% plot point estimate point-dropping averages

% point estimate selection
type = 'median';

% figure initialization
figure(...
    figopt,...
    'name',sprintf('point estimates (%s)',type),...
    'numbertitle','off');

% axes initialization
axes(...
    axesopt.default,...
    'xlim',[0,t_set(end)],...
    'xtick',sort([0;t_set]),...
    'ylim',[0,t_set(end)],...
    'ytick',sort([0;t_set]));
xlabel('Real time since T_2 onset (ms)');
ylabel('Decoded time since T_2 onset (ms)');

% plot unity line
plot(xlim,ylim,'--k');

% iterate through intensities
for ii = 1 : n_i
    i2_flags = concat_contrasts(nbdopt.test.trial_idcs) == i_set(ii);
    t2_flags = ismember(concat_stimuli(nbdopt.test.trial_idcs),t_set(stim2test_idcs(1:end)));
    concat_flags = ...
        t2_flags & ...
        i2_flags;
    
    pthat_avg = pthat_avgfuns{ff}(pthat.(type)(:,concat_flags));
    pthat_err = pthat_errfuns{ff}(pthat.(type)(:,concat_flags));
    
    t2bin_idcs = (1 : t_set(stim2test_idcs(end))) * psthbin;
    
    xpatch = [nbdopt.time(t2bin_idcs),fliplr(nbdopt.time(t2bin_idcs))];
    ypatch = [pthat_err(t2bin_idcs,1);flipud(pthat_err(t2bin_idcs,2))];
    patch(xpatch,ypatch,contrast_clrs(ii,:),...
        'facecolor',contrast_clrs(ii,:),...
        'edgecolor','none',...
        'facealpha',.15+(i_set(ii) == i_set(end))*.15);
    plot(nbdopt.time(t2bin_idcs),pthat_avg(t2bin_idcs),...
        'color',contrast_clrs(ii,:),...
        'linewidth',1.5);
    
    onset_idx = 1;
    plot(nbdopt.time(onset_idx),pthat_avg(onset_idx),...
        'color',contrast_clrs(ii,:),...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',8.5,...
        'markerfacecolor','w',...
        'markeredgecolor',contrast_clrs(ii,:));
    
    % iterate through stimuli
    for jj = 1 : stim2test_n
        offset_idx = find(nbdopt.time >= t_set(stim2test_idcs(jj)),1) - 1;
        plot(nbdopt.time(offset_idx),pthat_avg(offset_idx),...
            'color',contrast_clrs(ii,:),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',8.5,...
            'markerfacecolor',contrast_clrs(ii,:),...
            'markeredgecolor','k');
    end
end