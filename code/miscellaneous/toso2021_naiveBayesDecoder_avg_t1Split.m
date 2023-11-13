%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% training & test set conditions
for tt = 1 : n_t
    
    % training set conditions
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set(tt),[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
    
    % test set conditions
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),t_set(tt),[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
    
    n_conditions = conditions.test.n + conditions.train.n;
    
    % print training & test conditions
    fprintf('\nTRAINING CONDITIONS:\n');
    conditions.train.values
    fprintf('\nTEST CONDITIONS:\n');
    conditions.test.values
    
    %% run settings
    n_runs = 30;
    
    %% concatenation settings
    n_concatspercond = 2^7; % 2^8;
    n_concats = n_concatspercond * n_conditions;
    
    %% time settings
    roi = [-500,t_set(end)];
    roi_n_bins = range(roi) / psthbin;
    roi_time = linspace(roi(1),roi(2),roi_n_bins);
    
    %% construct spike rate tensor (time X neurons X concatenations)
    
    % data clearance
    clear concat_tensor P_tR;
    
    % preallocation
    P_tR = nan(roi_n_bins,roi_n_bins,conditions.test.n,n_runs);
    % concat_stimuli = nan(n_concats,n_runs);
    % concat_contrasts = nan(n_concats,n_runs);
    % concat_choices = nan(n_concats,n_runs);
    % concat_evalset = categorical(nan(n_concats,n_runs),[0,1],{'train','test'});
    
    % data type selection
    spike_data_field = 'FR';
    
    % iterate through runs
    for rr = 1 : n_runs
        
        % preallocation
        R = nan(roi_n_bins,n_neurons,n_conditions);
        %     concat_tensor = nan(roi_n_bins,n_neurons,n_concats);
        %     concat_stimuli = nan(n_concats,1);
        %     concat_contrasts = nan(n_concats,1);
        %     concat_choices = nan(n_concats,1);
        %     concat_evalset = categorical(nan(n_concats,1),[0,1],{'train','test'});
        
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
                %             concat_idcs = (1 : n_concatspercond) + ...
                %                 n_concatspercond * (kk - 1);
                %             concat_tensor(:,nn,concat_idcs) = aligned_spkrates(rand_idcs,:)';
                %             concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
                %             concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
                %             concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
                %             concat_evalset(concat_idcs) = 'train';
                
                r = aligned_spkrates(rand_idcs,:);
                r_mu = nanmean(r);
                nan_flags = isnan(r_mu);
                
                %             t_mat = roi_time' + (0 : n_concatspercond - 1);
                %             r_mat = r';
                %             nan_flags = isnan(r_mat);
                %             r_vec = r_mat(~nan_flags);
                %             t_vec = t_mat(~nan_flags);
                %             mdl = fit(t_vec,r_vec,'poly9',...
                %                 'robust','lar');
                %             mdl_poly = fit(t_vec,r_vec,'poly9');
                %             r_poly = max(0,mdl_poly(roi_time));
                
%                 mdl_poly = fit(roi_time(~nan_flags)',r_mu(~nan_flags)','poly9');
%                 r_poly = max(realmin,mdl_poly(roi_time));
%                 r_poly(nan_flags) = nan;
                
                mdl_spline = fit(...
                    roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
                    'smoothingparam',1e-6);
                r_spline = max(realmin,mdl_spline(roi_time));
                r_spline(nan_flags) = nan;
                
                %             p = polyfit(t_vec,r_vec,7);
                %             r_poly = max(0,polyval(p,roi_time));
                
                %             mdl = fit(roi_time(~nan_flags)',r_mu(~nan_flags)','poly9',...
                %                 'weights',sum(~isnan(r(:,~nan_flags))));
                %             r_polyw = max(0,mdl(roi_time));
                
%                 figure('position',[119.4000 53.8000 560 712.8000]);
%                 subplot(3,1,[1,2]);
%                 imagesc(r);
%                 subplot(3,1,3); hold on;
%                 plot(roi_time,r_mu);
%                 plot(roi_time,r_poly);
%                 plot(roi_time,r_spline);
                
                R(:,nn,kk) = r_spline;
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
                %             concat_idcs = (1 : n_concatspercond) + ...
                %                 n_concatspercond * (kk + conditions.train.n - 1);
                %             concat_tensor(:,nn,concat_idcs) = aligned_spkrates(rand_idcs,:)';
                %             concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
                %             concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
                %             concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
                %             concat_evalset(concat_idcs) = 'test';
                
                r = aligned_spkrates(rand_idcs,:);
                r_mu = nanmean(r);
                nan_flags = isnan(r_mu);
                
                %             t_mat = roi_time' + (0 : n_concatspercond - 1);
                %             r_mat = r';
                %             nan_flags = isnan(r_mat);
                %             r_vec = r_mat(~nan_flags);
                %             t_vec = t_mat(~nan_flags);
                %             mdl = fit(t_vec,r_vec,'poly9',...
                %                 'robust','lar');
                %             mdl_poly = fit(t_vec,r_vec,'poly9');
                %             r_poly = max(0,mdl_poly(roi_time));

%                 mdl_poly = fit(roi_time(~nan_flags)',r_mu(~nan_flags)','poly9');
%                 r_poly = max(realmin,mdl_poly(roi_time));
%                 r_poly(nan_flags) = nan;
                
                mdl_spline = fit(...
                    roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
                    'smoothingparam',1e-6);
                r_spline = max(realmin,mdl_spline(roi_time));
                r_spline(nan_flags) = nan;
                
                %             p = polyfit(roi_time,r_mu,7);
                %             r_poly = max(0,polyval(p,roi_time));
                
                %             mdl = fit(roi_time(~nan_flags)',r_mu(~nan_flags)','poly9',...
                %                 'weights',sum(~isnan(r(:,~nan_flags))));
                %             r_polyw = max(0,mdl(roi_time));
                
%                 figure('position',[119.4000 53.8000 560 712.8000]);
%                 subplot(3,1,[1,2]);
%                 imagesc([roi],[],r);
%                 subplot(3,1,3); hold on;
%                 plot(roi_time,r_mu);
%                 plot(roi_time,r_poly);
%                 plot(roi_time,r_spline);
                
                R(:,nn,kk+conditions.train.n) = r_spline;
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
        P_tR(:,:,:,rr) = naivebayestimedecoder(R,nbdopt);
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
        'xlim',[roi(1),t_set(end)],...
        'xtick',unique([roi';0;t_set]),...
        'ylim',[roi(1),t_set(end)],...
        'ytick',unique([roi';0;t_set]));
    xlabel('Real time since S_2 onset (ms)');
    ylabel('Decoded time since S_2 onset (ms)');
    
    % posterior subtraction
    p_contrast_min = avgfun(P_tR(:,:,1,:),4);
    p_contrast_max = avgfun(P_tR(:,:,3,:),4);
    p_diff = p_contrast_max - p_contrast_min;
    imagesc(xlim,ylim,p_diff',...
        [-1,1] * n_t / n_tbins * 5);
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
        xlabel(sps(ii),'Real time since S_2 onset (ms)');
        ylabel(sps(ii),'Decoded time since S_2 onset (ms)');
    end
    set(sps,...
        axesopt.default,...
        'xlim',[roi(1),t_set(end)],...
        'xtick',unique([roi';0;t_set]),...
        'ylim',[roi(1),t_set(end)],...
        'ytick',unique([roi';0;t_set]));
    linkaxes(sps);
    
    clims = quantile(P_tR,[0,.999],'all')';
    % iterate through contrast conditions
    for ii = 1 : n_contrasts
        title(sps(ii),sprintf('%s = %.0f mm/s',...
            contrast_lbl,contrast_set(ii)),...
            'fontsize',10);
        p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
        p_cond(isnan(p_cond)) = max(clims);
        imagesc(sps(ii),xlim,ylim,p_cond',clims);
        plot(sps(ii),xlim,ylim,'--w');
    end
    
    %% plot superimposed contrast-split posterior averages
    figure(...
        'name','posterior averages',...
        'numbertitle','off');
    axes(...
        axesopt.default,...
        'xlim',[roi(1),t_set(end)],...
        'xtick',unique([roi';0;t_set]),...
        'ylim',[roi(1),t_set(end)],...
        'ytick',unique([roi';0;t_set]));
    xlabel('Real time since S_2 onset (ms)');
    ylabel('Decoded time since S_2 onset (ms)');
    
    % color limits
    clims = quantile(avgfun(P_tR,4),[0,1],'all')';
    
    % iterate through contrast conditions
    for ii = 1 : n_contrasts
        p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
        p_cond(p_cond < clims(1)) = clims(1);
        p_cond(p_cond > clims(2)) = clims(2);
        p_cond(isnan(p_cond)) = max(p_cond(:));
        p_patch = mat2patch(p_cond,xlim,ylim);
        patch(p_patch,...
            'facevertexalphadata',p_patch.facevertexcdata,...
            'edgecolor','none',...
            'facecolor',contrast_clrs(ii,:),...
            'facealpha','flat');
    end
    
    % plot identity line
    plot(xlim,ylim,'--k');
end