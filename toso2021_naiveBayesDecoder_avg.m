%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% training & test set conditions
t_idx = 1;

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
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
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
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
end

n_conditions = conditions.test.n + conditions.train.n;

% print training & test conditions
fprintf('\nTRAINING CONDITIONS:\n');
conditions.train.values
fprintf('\nTEST CONDITIONS:\n');
conditions.test.values

%% run settings
n_runs = 10;

%% concatenation settings
n_concatspercond = 2^8; % 2^8;
n_concats = n_concatspercond * n_conditions;

%% time settings
roi = [-500,t_set(end)]; % [-500,t_set(end)] !!!!!!!!!!!
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
            
%             time_mat = roi_time' + (0 : n_concatspercond - 1);% * range(roi);
%             nan_flags = isnan(r);
%             r_mat = r(~nan_flags)';
%             time_mat = time_mat(~nan_flags)';
%             p = polyfit(time_mat(:),r_mat(:),10);
%             r_poly = polyval(p,roi_time);
            
            p = polyfit(roi_time,r_mu,10);
            r_poly = max(0,polyval(p,roi_time));
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);

            R(:,nn,kk) = r_poly;
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
            
            p = polyfit(roi_time,r_mu,10);
            r_poly = max(0,polyval(p,roi_time));

            R(:,nn,kk+conditions.train.n) = r_poly;
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
    'xlim',roi,...
    'xtick',sort([0;t_set]),...
    'ylim',roi,...
    'ytick',sort([0;t_set]));
xlabel('Real time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% posterior subtraction
p_contrast_min = avgfun(P_tR(:,:,1,:),4);
p_contrast_max = avgfun(P_tR(:,:,3,:),4);
p_diff = p_contrast_max - p_contrast_min;
imagesc(roi,roi,p_diff',...
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
end
set(sps,...
    axesopt.default,...
    'xlim',roi,...
    'xtick',[],...sort([0;t_set]),...
    'ylim',roi,...
    'ytick',[]); % sort([0;t_set]));
linkaxes(sps);

clims = quantile(P_tR,[0,.999],'all')';
% iterate through contrast conditions
for ii = 1 : n_contrasts
    title(sps(ii),sprintf('%s = %.0f mm/s',...
        contrast_lbl,contrast_set(ii)),...
        'fontsize',10);
    p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
    imagesc(sps(ii),roi,roi,p_cond',clims);
    plot(sps(ii),xlim,ylim,'--w');
end