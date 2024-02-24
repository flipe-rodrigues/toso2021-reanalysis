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
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),t_set(t2_mode_idx),[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),[],t_set([1:t2_mode_idx+1,end]),...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
%     conditions.train = intersectconditions(...
%         't1',t1(valid_flags),[],[],...
%         'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
%         't2',t2(valid_flags),[],t_set([1:t2_mode_idx+1,end]),...
%         'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
%         'choice',choice(valid_flags),[],[]);
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
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
        't1',t1(valid_flags),t_set([2,end-1]),[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),t_set(t2_mode_idx),[],...
        'i2',i2(valid_flags),i_set(i1_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set,[],...
        't2',t2(valid_flags),t_set(t2_mode_idx+1),[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
%     conditions.test = intersectconditions(...
%         't1',t1(valid_flags),[],[],...
%         'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
%         't2',t2(valid_flags),t_set(t2_mode_idx+1),[],...
%         'i2',i2(valid_flags),i_set,[],...
%         'choice',choice(valid_flags),[],[]);
	conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
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
n_runs = 10;

%% concatenation settings
n_concatspercond = 2^6;
n_concats = n_concatspercond * n_conditions;

%% time settings
roi = [-500,max(conditions.train.values.t2)];
roi_n_bins = range(roi) / psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

%% selection criteria
n_trial_cutoff = 4;

%% construct spike rate tensor (time X neurons X concatenations)

% preallocation
scale_factors = nan(n_neurons,conditions.test.n,n_runs);

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
        xval_train_trials = cell(conditions.train.n,conditions.test.n);
        
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
            if n_flagged_trials <= n_trial_cutoff
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_rates = data.SDF(trial_flags,:)';
            
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
                    xval_train_trials{kk,ii} = xval_trials(xval_train_idcs);
                end
                
                % concatenate sub-sampled training sets across test conditions
                train_idcs = find(ismember(...
                    flagged_trials,vertcat(xval_train_trials{kk,:})));
            else
                
                % train using all trials in the remaining conditions
                train_idcs = 1 : n_flagged_trials;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(train_idcs,n_concatspercond,true);

            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
%             nan_flags = isnan(r_mu);

%             mdl_spline = fit(...
%                 roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
%                 'smoothingparam',1e-6);
%             r_spline = max(1e-3,mdl_spline(roi_time));
%             r_spline(nan_flags) = nan;
            
%             t_mat = repmat(roi_time,n_concatspercond,1)';
%             r_mat = r';
%             r_vec = r_mat(~isnan(r_mat));
%             t_vec = t_mat(~isnan(r_mat));
%             mdl_poly = fit(t_vec,r_vec,'poly9');
%             r_poly = max(1e-3,mdl_poly(roi_time));
%             r_poly(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
%             plot(roi_time,r_spline);
%             plot(roi_time,(r_poly+r_spline)/2);
            
            R(:,nn,kk) = r_mu; % nanmean(aligned_spkrates);
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
            if n_flagged_trials <= n_trial_cutoff
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_rates = data.SDF(trial_flags,:)';
            
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
            aligned_spkrates = spike_rates;
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [roi_n_bins,n_flagged_trials])';
            
            % handle cross-validation
            test_idcs = find(~ismember(...
                flagged_trials,vertcat(xval_train_trials{:,kk})));
            if isempty(test_idcs)
                continue;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(test_idcs,n_concatspercond,true);

            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
%             nan_flags = isnan(r_mu);
            
%             mdl_spline = fit(...
%                 roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
%                 'smoothingparam',1e-6);
%             r_spline = max(1e-3,mdl_spline(roi_time));
%             r_spline(nan_flags) = nan;
            
%             t_mat = repmat(roi_time,n_concatspercond,1)';
%             r_mat = r';
%             r_vec = r_mat(~isnan(r_mat));
%             t_vec = t_mat(~isnan(r_mat));
%             mdl_poly = fit(t_vec,r_vec,'poly9');
%             r_poly = max(1e-3,mdl_poly(roi_time));
%             r_poly(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
%             plot(roi_time,r_spline);
            
            R(:,nn,kk+conditions.train.n) = r_mu; % nanmean(aligned_spkrates);
        end
    end
    
    %% compute scaling factors
    
    % temporal warping settings
    warpopt.bounds = [.25,2]; % 1 + ([-25,50]+[-1,1]*75*.1) / 100;
    warpopt.evalbounds = [0,inf];
    warpopt.upsfactor = 1; % ceil(3 * range(warpopt.bounds));
    warpopt.lambda = .25;
    warpopt.clrs = contrast_clrs;
    warpopt.verbose = true;
    warpopt.center = 0;
    warpopt.normalize = 0;
    warpopt.pullapart = 0;
    warpopt.plot = 0;
    
    % compute observation weights
    time_mat = repmat(roi_time,n_total_trials,1);
    weights = sum(time_mat <= t2);
    weights = weights ./ max(weights);
    warpopt.weights = weights .^ 1;
    
    % fit temporal scaling factors
    tic
    [~,scale_factors(:,:,rr)] = fittemporalscalingfactor_temporalscaleonly(...
        R(:,:,1),R(:,:,2:end),roi_time,warpopt);
    toc
end

%% plot response dilation

% set dilation boundaries
dilation_bounds = [-1,1] * 50;

% figure initialization
fig = figure(figopt,...
    'name',sprintf('response_dilation_%s',contrast_lbl));
axes(axesopt.default,...
    'xlim',[1,conditions.test.n] + [-1,1]*.65,...
    'xtick',1:conditions.test.n,...
    'xticklabel',num2cell(conditions.test.values.(contrast_str)),...
    'ylim',dilation_bounds+[-1,1]*range(dilation_bounds)*.05,...
    'ytick',linspace(dilation_bounds(1),dilation_bounds(2),5),...
    'layer','bottom',...
    'plotboxaspectratio',[conditions.test.n+1,8,1]);
xlabel(sprintf('%s (%s)',contrast_lbl,contrast_units));
ylabel('Response dilation (%)');

% convert scaling factors to dilations
dilations = (nanmean(scale_factors,3) - 1) * 100;

% zero line
plot(xlim,[1,1]*0,...
    'color','k',...
    'linestyle',':');

% graphical object preallocation
p = gobjects(3,conditions.test.n);

% bar settings
barwidth = .75;

% dilation selection settings
lowerbound = min(ylim);
upperbound = max(ylim);
% lowerbound = -inf;
% upperbound = +inf;

% average function selection
avgfun = @(x,d) nanmedian(x,d);
errfun = @(x) diff(quantile(x,3));
% avgfun = @(x) nanmean(x);
% errfun = @(x) [1,1] * nanstd(x) / sqrt(numel(x));

% compute the reference average
fov_flags = all(...
    dilations >= lowerbound & ...
    dilations <= upperbound,2);
offset = avgfun(dilations(fov_flags,:),[1,2]);

% 
bound_flags = all(...
    dilations ~= (warpopt.bounds(1) - 1) * 100 & ...
    dilations ~= (warpopt.bounds(2) - 1) * 100,2);

%
dilation_flags = ...
    ...fov_flags & ...
    ...bound_flags & ...
    true(n_neurons,1);
n_flaggeddilations = sum(dilation_flags)
dilation_edges = linspace(...
    dilation_bounds(1)*1.25,dilation_bounds(2)*1.25,n_flaggeddilations);

% iterate through conditions
for kk = 1 : conditions.test.n

    % mean & standard error across animals
    cond_avg = avgfun(dilations(dilation_flags,kk),1) - offset;
    cond_err = errfun(dilations(dilation_flags,kk));
    
    % plot summary of threshold dilation
    xpatch = kk + [-1,1,1,-1] * barwidth / 2;
    ypatch = cond_avg + [[-1,-1] * cond_err(1),[+1,+1] * cond_err(2)];

    p(1,kk) = plot(kk+[-1,1]*.65*barwidth,[1,1]*cond_avg,...
        'color',contrast_clrs(conditions.test.idcs.(contrast_str)(kk),:),...
        'linestyle','-',...
        'linewidth',2.5);
    p(2,kk) = plot(kk+[-1,1]*.5*barwidth,[1,1]*cond_avg,...
        'color','w',...
        'linestyle','-',...
        'linewidth',5);    
    p(3,kk) = patch(xpatch,ypatch,contrast_clrs(conditions.test.idcs.(contrast_str)(kk),:),...
        'edgecolor','w',...
        'facealpha',1,...
        'linewidth',1);
    
    % plot ground truth
    plot(xlim,[1,1]*(contrast_set(conditions.test.idcs.(contrast_str)(kk))-1)*100,...
        'color',contrast_clrs(conditions.test.idcs.(contrast_str)(kk),:),...
        'linestyle',':',...
        'linewidth',1);
end

% iterate through conditions
for kk = 1 : conditions.test.n
    [dilation_pdf,~] = ksdensity(...
        dilations(dilation_flags,kk) - offset,dilation_edges);
    dilation_pdf = .2 * ...
        (dilation_pdf - min(dilation_pdf)) ./ range(dilation_pdf);
    noise = randn(n_flaggeddilations,1) .* dilation_pdf';

    % plot animal threshold dilation
    grapeplot(kk+noise,...
        sort(dilations(dilation_flags,kk)) - offset,...
        'marker','o',...
        'markersize',5,...
        'markeredgecolor',contrast_clrs(conditions.test.idcs.(contrast_str)(kk),:),...
        'markerfacecolor',[1,1,1],...
        'linewidth',1);
end

% ui restacking
sp = p(1:2,:); 
uistack(sp(:),'top');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% compute response stretch

% preallocation
stretch = nan(n_neurons,1);

% design matrix
x = conditions.test.values.(contrast_str);% - contrast_set(contrast_mode_idx);
X = [ones(conditions.test.n,1),x];

% iterate through units
for uu = 1 : n_neurons

    % linear regression
    thetas = X \ (dilations(uu,:)' - offset);
    stretch(uu) = thetas(end);
end

% one-sample t-test
[~,pvalue,~,stats.ttest] = ttest(stretch);
stats.ttest.pvalue = pvalue;

%% plot response stretch
stretch_bounds = [-1,1] * 1.5;
count_bounds = [0,60];

% figure initialization
fig = figure(figopt,...
    'name',sprintf('response_stretch_%s',contrast_lbl));
axes(...
    axesopt.default,...
    'plotboxaspectratio',[1,3,1],...
    'xlim',count_bounds,...
    'xtick',count_bounds,...
    'xdir','reverse',...
    'ylim',stretch_bounds+[-1,1]*.05*range(stretch_bounds),...
    'ytick',linspace(stretch_bounds(1),stretch_bounds(2),5),...
    'yaxislocation','right');
xlabel('Count');
ylabel(sprintf('Response stretch (%% / %s)',contrast_units),...
    'rotation',-90,...
    'verticalalignment','bottom');

% zero line
plot(xlim,[1,1]*0,...
    'color','k',...
    'linestyle',':');

% distribution settings
nbins = 41;
bin_edges = linspace(min(yticks),max(yticks),nbins);
bin_counts = histcounts(stretch,bin_edges);

% plot stretch distribution
histogram(...
    'binedges',bin_edges,...
    'bincounts',bin_counts,...
    'orientation','horizontal',...
    'edgecolor','k',...
    'facecolor','k',...
    'facealpha',.5,...
    'linestyle','none');
stairs([0,bin_counts],bin_edges,...
    'color','k',...
    'linewidth',1.5);

% compute mean stretch
fov_flags = ...
    stretch >= min(ylim) & ...
    stretch <= max(ylim);
stretch_avg = nanmean(stretch(fov_flags));

% assess population significance
if pvalue <= .01
    str = '**';
elseif pvalue <= .05
    str = '*';
else
    str = '^{n.s.}';
end
text(max(xlim)*.9,.0375*range(ylim),str,...
    'color','k',...
    'fontname',axesopt.default.fontname,....
    'fontsize',axesopt.default.fontsize*1.25,...
    'horizontalalignment','center');

% test affordance
p = plot([1,1]*max(xlim)*.9,[0,1]*stretch_avg,...
    'color','k',...
    'linewidth',1.5);
uistack(p,'bottom');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end