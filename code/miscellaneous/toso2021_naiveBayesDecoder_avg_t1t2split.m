%% initialization
if ~exist('data','var')
    toso2021_wrapper;
    close all;
end

%% run settings
n_runs = 50;

%% concatenation settings
n_concats_max = 2^7;

%% time settings
roi = [-500,t_set(end)];
roi2plot = [-240,t_set(t2_mode_idx+1)];
roi_n_bins = range(roi) / psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

%% slice settings
slice_times = unique([roi2plot(1)/2;0;t_set(1:end-2)]);
n_slices = numel(slice_times);
slice_offsets = (n_slices - 1 : -1 : 0) * .05 + .015;
slice_clrs = gray(n_slices);

%% T1- and T2-specific decoding

% preallocation
P_tR_avgs = nan(roi_n_bins,roi_n_bins,n_contrasts,n_s_pairs);
map_avgs = nan(roi_n_bins,n_contrasts,n_s_pairs);

% iterate through S1-S2 pairs
for ss = 1 : n_s_pairs
    
    % training set conditions
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),s_pairset(ss,1),[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),s_pairset(ss,2),[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
    
    % test set conditions
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),s_pairset(ss,1),[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),s_pairset(ss,2),[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
    n_conditions = conditions.test.n + conditions.train.n;
    
    % print training & test conditions
    fprintf('\nTRAINING CONDITIONS:\n');
    conditions.train.values
    fprintf('\nTEST CONDITIONS:\n');
    conditions.test.values
    
    %% construct spike rate tensor (time X neurons X concatenations)
    
    % data clearance
    clear concat_tensor P_tR;
    
    % preallocation
    P_tR = nan(roi_n_bins,roi_n_bins,conditions.test.n,n_runs);
    map = nan(roi_n_bins,conditions.test.n,n_runs);
    
    % temporal smoothing kernel
    gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);
    
    % iterate through runs
    for rr = 1 : n_runs
        
        % preallocation
        R = nan(roi_n_bins,n_neurons,n_conditions);
        
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
                        n_train_trials = round(n_xval_trials * 2 / 3);
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
                rand_idcs = randsample(train_idcs,n_concats_max,true);
                
                r = aligned_spkrates(rand_idcs,:);
                r_mu = nanmean(r);
                nan_flags = isnan(r_mu);
                
                r_gauss = nanconv2(r_mu,1,gauss_kernel.pdf);
                r_gauss(nan_flags) = nan;
                
                R(:,nn,kk) = r_gauss;
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
                rand_idcs = randsample(test_idcs,n_concats_max,true);
                
                r = aligned_spkrates(rand_idcs,:);
                r_mu = nanmean(r);
                nan_flags = isnan(r_mu);
                
                r_gauss = nanconv2(r_mu,1,gauss_kernel.pdf);
                r_gauss(nan_flags) = nan;
                
                R(:,nn,kk+conditions.train.n) = r_gauss;
            end
        end
        
        %% naive bayes decoder
        nbdopt = struct();
        nbdopt.n_xpoints = 100;
        nbdopt.time = roi_time;
        nbdopt.train.trial_idcs = 1 : conditions.train.n;
        nbdopt.train.n_trials = numel(nbdopt.train.trial_idcs);
        nbdopt.test.trial_idcs = (1 : conditions.test.n) + conditions.train.n;
        nbdopt.test.n_trials = numel(nbdopt.test.trial_idcs);
        nbdopt.assumepoissonmdl = true;
        nbdopt.verbose = true;
        
        tic
        [P_tR(:,:,:,rr),P_Rt,pthat,neurons] = naivebayestimedecoder(R,nbdopt);
        map(:,:,rr) = pthat.mode;
        toc
    end

    %% average across runs
    P_tR_avgs(:,:,:,ss) = nanmean(P_tR,4);
    map_avgs(:,:,ss) = nanmean(map,3);
end

%% compute stimulus pair weights

% preallocation
s_pair_counts = nan(roi_n_bins,roi_n_bins,n_s_pairs);

% iterate through stimulus pairs
for ss = 1 : n_s_pairs
    s_pair_flags = ...
        t1 == s_pairset(ss,1) & ...
        t2 == s_pairset(ss,2);
    t2_flags = roi_time <= s_pairset(ss,2);
    s_pair_counts(t2_flags,t2_flags,ss) = sum(s_pair_flags);
end

% normalization
s_pair_weights = s_pair_counts ./ nansum(s_pair_counts,3);

% nan correction
s_pair_weights(isnan(s_pair_weights)) = 0;

%% compute average across stimulus pair

% preallocation
P_tR_avg = zeros(roi_n_bins,roi_n_bins,n_contrasts);
map_avg = zeros(roi_n_bins,n_contrasts);

% iterate through stimulus pairs
for ss = 1 : n_s_pairs
    
    % update weighted averages
    P_tR_ss = P_tR_avgs(:,:,:,ss);
%     P_tR_ss = P_tR_ss ./ nansum(P_tR_ss,2);
    P_tR_ss(isnan(P_tR_ss)) = 0;
    P_tR_avg = P_tR_avg + s_pair_weights(:,:,ss) .* P_tR_ss;
    map_ss = map_avgs(:,:,ss);
    map_ss(isnan(map_ss)) = 0;
    map_avg = map_avg + s_pair_weights(:,1,ss) .* map_ss;
    
    % figure initialization
    figure(...
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
    title(sprintf('T1=%i, T2=%i',s_pairset(ss,1),s_pairset(ss,2)));
    xlabel('Time since S_2 onset (ms)');
    ylabel('Decoded time since S_2 onset (ms)');
    
    % convert from tensor to rgb
%     P_tR_ss = P_tR_ss ./ nansum(P_tR_ss,2);
%     P_tR_ss = min(P_tR_ss,quantile(P_tR_ss,.995,'all'));
    P = tensor2rgb(permute(P_tR_ss,[2,1,3]),contrast_clrs);
    imagesc(roi,roi,P);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
    
        % plot MAP
        plot(roi_time(roi_time<=s_pairset(ss,2)),...
            map_ss(roi_time<=s_pairset(ss,2),ii),...
            'color',contrast_clrs(ii,:),...
            'linewidth',.5);
    end
    
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
    p_contrast_min = P_tR_ss(:,:,1);
    p_contrast_max = P_tR_ss(:,:,end);
    p_contrast_min = p_contrast_min - min(p_contrast_min,[],[1,2]);
    p_contrast_max = p_contrast_max - min(p_contrast_max,[],[1,2]);
    p_diff = p_contrast_max - p_contrast_min;
    imagesc(roi,roi,p_diff',[-1,1] * n_t / n_tbins * 5);
    
    % zero lines
    plot(xlim,[1,1]*0,':k');
    plot([1,1]*0,ylim,':k');
end

% normalization
% P_tR_avg = P_tR_avg ./ nansum(P_tR_avg,2);

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
% P_tR_avg = squeeze(avgfun(P_tR_avgs,4));
% P_tR_avg(roi_time < roi2plot(1) | roi_time > roi2plot(2),:,:) = nan;
% P_tR_avg(:,roi_time < roi2plot(1) | roi_time > roi2plot(2),:) = nan;
% P_tR_avg = P_tR_avg ./ nansum(P_tR_avg,2);
% P_tR_avg(find(roi_time>0,1):find(roi_time>10,1),:,3) = max(P_tR_avg,[],'all');
% P_tR_avg = min(P_tR_avg,quantile(P_tR_avg,.9995,'all'));
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
p_contrast_min = P_tR_avg(:,:,1);
p_contrast_max = P_tR_avg(:,:,end);
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
    
%     % should recompute the posterior means & stds....
%     [~,map_idcs] = max(P_tR_avg(:,:,ii),[],2);
%     p_map = roi_time(map_idcs);
%     
%     % posterior median
%     med_flags = [false(roi_n_bins,1),...
%         diff(cumsum(P_tR_avg(:,:,ii),2) > .5,1,2) == 1];
%     [~,med_idcs] = max(med_flags,[],2);
%     p_med = roi_time(med_idcs);
%     
%     p = P_tR_avg(:,:,ii);
%     p(isnan(p)) = 0;
%     p_mu = roi_time * p';
% 
%     plot(roi_time,p_map,...
%         'color',contrast_clrs(ii,:),...
%         'linewidth',1.5);
%     continue;
    
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