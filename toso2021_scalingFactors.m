%% initialization
if ~exist('data','var')
    toso2021_wrapper;
    close all;
end

%% training & test set conditions

% training set conditions
if strcmpi(contrast_str,'t1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set(t2_mode_idx+2),[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'choice')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],t_set([1,end]),...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),1,[]);
elseif strcmpi(contrast_str,'correct')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'correct',correct(valid_flags),1,[]);
end

% test set conditions
if strcmpi(contrast_str,'t1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set,[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set(t2_mode_idx),[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'choice')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[0,1],[]);
elseif strcmpi(contrast_str,'correct')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'correct',correct(valid_flags),[0,1],[]);
end

n_conditions = conditions.test.n + conditions.train.n;

% print training & test conditions
fprintf('\nTRAINING CONDITIONS:\n');
conditions.train.values
fprintf('\nTEST CONDITIONS:\n');
conditions.test.values

%% run settings
n_runs = 1;

%% concatenation settings
n_concatspercond = 2^6;
n_concats = n_concatspercond * n_conditions;

%% time settings
roi = [-500,max(conditions.train.values.t2)];
roi_n_bins = range(roi) / psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

%% construct spike rate tensor (time X neurons X concatenations)

% data clearance
clear concat_tensor P_tR;

% preallocation
P_tR = nan(roi_n_bins,roi_n_bins,conditions.test.n,n_runs);
MAP = nan(roi_n_bins,conditions.test.n,n_runs);

% data type selection
spike_data_field = 'FR';

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
            rand_idcs = randsample(train_idcs,n_concatspercond,true);

            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
            nan_flags = isnan(r_mu);

%             mdl_spline = fit(...
%                 roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
%                 'smoothingparam',1e-6);
%             r_spline = max(1e-3,mdl_spline(roi_time));
%             r_spline(nan_flags) = nan;
            
            t_mat = repmat(roi_time,n_concatspercond,1)';
            r_mat = r';
            r_vec = r_mat(~isnan(r_mat));
            t_vec = t_mat(~isnan(r_mat));
            mdl_poly = fit(t_vec,r_vec,'poly9');
            r_poly = max(1e-3,mdl_poly(roi_time));
            r_poly(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
%             plot(roi_time,r_spline);
%             plot(roi_time,(r_poly+r_spline)/2);
            
            R(:,nn,kk) = r_mu;
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

            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
            nan_flags = isnan(r_mu);
            
%             mdl_spline = fit(...
%                 roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
%                 'smoothingparam',1e-6);
%             r_spline = max(1e-3,mdl_spline(roi_time));
%             r_spline(nan_flags) = nan;
            
            t_mat = repmat(roi_time,n_concatspercond,1)';
            r_mat = r';
            r_vec = r_mat(~isnan(r_mat));
            t_vec = t_mat(~isnan(r_mat));
            mdl_poly = fit(t_vec,r_vec,'poly9');
            r_poly = max(1e-3,mdl_poly(roi_time));
            r_poly(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
%             plot(roi_time,r_spline);
            
            R(:,nn,kk+conditions.train.n) = r_mu;
        end
    end
    
    %% compute scaling factors
    
    % temporal warping settings
    warpopt.bounds = [.5,1.5]; % 1 + ([-25,50]+[-1,1]*75*.1) / 100;
    warpopt.evalbounds = [0,inf];
    warpopt.upsfactor = ceil(3 * range(warpopt.bounds));
    warpopt.lambda = .0;
    warpopt.clrs = contrast_clrs;
    warpopt.verbose = true;
    warpopt.center = 0;
    warpopt.normalize = 0;
    warpopt.plot = 0;
    warpopt.weights = roi_time .^ 0;
    tic
    [~,scale_factors] = fittemporalscalingfactor(...
        R(:,:,1),R(:,:,2:end),roi_time,warpopt);
    toc
end

%% plot response dilation

% set dilation boundaries
dilation_bounds = [-1,1] * 50;

% figure initialization
figure(figopt.default,...
    'name',sprintf('response dilation (%s)',contrast_lbl));
axes(axesopt.default,...
    'xlim',[1,n_contrasts] + [-1,1]*.65,...
    'xtick',1:n_contrasts,...
    'xticklabel',num2cell(contrast_set),...
    'ylim',dilation_bounds+[-1,1]*range(dilation_bounds)*.05,...
    'ytick',linspace(dilation_bounds(1),dilation_bounds(2),5),...
    'plotboxaspectratio',[n_contrasts,8,1]);
xlabel(sprintf('%s (%s)',contrast_lbl,contrast_units));
ylabel('Response dilation (%)');

% convert scaling factors to dilations
dilations = (scale_factors - 1) * 100;

% zero line
plot(xlim,[1,1]*0,...
    'color','k',...
    'linestyle',':');

% graphical object preallocation
p = gobjects(3,n_contrasts);

% bar settings
barwidth = .75;

% dilation selection settings
lowerbound = min(ylim);
upperbound = max(ylim);
lowerbound = -inf;
upperbound = +inf;

% average function selection
avgfun = @(x,d) nanmedian(x,d);
errfun = @(x) diff(quantile(x,3));
% avgfun = @(x) nanmean(x);
% errfun = @(x) [1,1] * nanstd(x) / sqrt(numel(x));

% compute the reference average
fov_flags = any(...
    dilations >= lowerbound & ...
    dilations <= upperbound,2);
n_fovunits = sum(fov_flags)
fov_edges = linspace(...
    dilation_bounds(1)*1.25,dilation_bounds(2)*1.25,n_fovunits);
offset = avgfun(dilations(fov_flags,:),[1,2]);

% iterate through conditions
for kk = 1 : n_contrasts

    % mean & standard error across animals
    cond_avg = avgfun(dilations(fov_flags,kk),1) - offset;
    cond_err = errfun(dilations(fov_flags,kk));
    
    % plot summary of threshold dilation
    xpatch = kk + [-1,1,1,-1] * barwidth / 2;
    ypatch = cond_avg + [[-1,-1] * cond_err(1),[+1,+1] * cond_err(2)];

    p(1,kk) = plot(kk+[-1,1]*.65*barwidth,[1,1]*cond_avg,...
        'color',contrast_clrs(kk,:),...
        'linestyle','-',...
        'linewidth',2.5);
    p(2,kk) = plot(kk+[-1,1]*.5*barwidth,[1,1]*cond_avg,...
        'color','w',...
        'linestyle','-',...
        'linewidth',5);    
    p(3,kk) = patch(xpatch,ypatch,contrast_clrs(kk,:),...
        'edgecolor','w',...
        'facealpha',1,...
        'linewidth',1);
    
    % plot ground truth
    plot(xlim,[1,1]*(contrast_set(kk)-1)*100,...
        'color',contrast_clrs(kk,:),...
        'linestyle',':',...
        'linewidth',1);
end

% iterate through conditions
for kk = 1 : n_contrasts
    [dilation_pdf,~] = ksdensity(...
        dilations(fov_flags,kk) - offset,fov_edges);
    dilation_pdf = .2 * ...
        (dilation_pdf - min(dilation_pdf)) ./ range(dilation_pdf);
    noise = randn(n_fovunits,1) .* dilation_pdf';

    % plot animal threshold dilation
    grapeplot(kk+noise,...
        sort(dilations(fov_flags,kk)) - offset,...
        'marker','o',...
        'markersize',5,...
        'markeredgecolor',contrast_clrs(kk,:),...
        'markerfacecolor','w',...
        'linewidth',1);
end

% ui restacking
sp = p(1:2,:); 
uistack(sp(:),'top');