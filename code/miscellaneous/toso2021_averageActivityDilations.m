%% initialization
if ~exist('data','var')
    toso2021_wrapper;
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
        't2',t2(valid_flags),[],t_set([1:t2_mode_idx+1,end]),...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
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
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set,[],...
        't2',t2(valid_flags),t_set(t2_mode_idx+1),[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
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

%% selection criteria
n_trial_cutoff = 0;

%% window settings
window_dur = t_set(t2_mode_idx);
t2_flags = t2 >= window_dur;

%% ROI settings

% initialization
roi_win = struct();
roi_lbl = struct();

% roi definitions
roi_win.pre = [-window_dur,0];
roi_win.post = [0,window_dur];

% roi labels
roi_lbl.pre = 'Pre-S2 onset';
roi_lbl.post = 'Post-S2 onset';

% roi length
roi_length.pre = range(roi_win.pre) / psthbin;
roi_length.post = range(roi_win.post) / psthbin;

% parse roi epochs
roi_epochs = fieldnames(roi_win);
n_epochs = numel(roi_epochs);

%% construct spike rate tensor (time X neurons X concatenations)

% preallocation
scale_factors = nan(n_neurons,conditions.test.n,n_runs,n_epochs);

% iterate through epochs
for ee = 1 : n_epochs
    epoch = roi_epochs{ee};
    
    % iterate through runs
    for rr = 1 : n_runs
        
        % preallocation
        R = nan(n_neurons,n_conditions);
        
        % iterate through units
        for nn = 1 : n_neurons
            progressreport(nn,n_neurons,...
                sprintf('fetching spikes (run %i/%i, %s)',rr,n_runs,epoch));
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
                    ...t2_flags & ...
                    condition_flags;
                flagged_trials = find(trial_flags);
                n_flagged_trials = numel(flagged_trials);
                if n_flagged_trials <= n_trial_cutoff
                    continue;
                end
                
                % fetch spike counts & compute spike rates
                spike_counts = data.FR(trial_flags,:);
                spike_rates = ...
                    conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
                
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        spike_rates = downsamplecounts(spike_counts,25);
        spike_rates = spike_rates(:,validtime_flags)';
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                % S2-onset-aligned spike rates
                alignment = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags) + ...
                    isi;
                alignment_flags = ...
                    valid_time >= alignment + roi_win.(epoch)(1) & ...
                    valid_time < alignment + t2(trial_flags);
                chunk_flags = ...
                    valid_time >= alignment + roi_win.(epoch)(1) & ...
                    valid_time < alignment + roi_win.(epoch)(2);
                aligned_spkrates = spike_rates;
                aligned_spkrates(~alignment_flags') = nan;
                aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                    [roi_length.(epoch),n_flagged_trials])';
                
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
                    for xx = xval_condition_idcs
                        
                        % detect overlap between training & test trials
                        feature_flags = false(n_flagged_trials,1);
                        for ff = 1 : conditions.test.features.n
                            feature_lbl = conditions.test.features.labels{ff};
                            feature = eval(feature_lbl);
                            feature_flags(:,ff) = any(ismember(...
                                feature(flagged_trials),...
                                conditions.test.values.(feature_lbl)(xx,:)),2);
                        end
                        xval_trial_flags = all(feature_flags,2);
                        
                        % split the conflicting trials into training & test subsets
                        xval_trials = flagged_trials(xval_trial_flags);
                        n_xval_trials = numel(xval_trials);
                        n_train_trials = round(n_xval_trials * 2 / 3);
                        xval_train_idcs = randperm(n_xval_trials,n_train_trials);
                        xval_train_trials{kk,xx} = xval_trials(xval_train_idcs);
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
                R(nn,kk) = nanmean(r,[1,2]);
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
                    ...t2_flags & ...
                    condition_flags;
                flagged_trials = find(trial_flags);
                n_flagged_trials = numel(flagged_trials);
                if n_flagged_trials <= n_trial_cutoff
                    continue;
                end
                
                % fetch spike counts & compute spike rates
                spike_counts = data.FR(trial_flags,:);
                spike_rates = ...
                    conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
                
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        spike_rates = downsamplecounts(spike_counts,25);
        spike_rates = spike_rates(:,validtime_flags)';
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
                % S2-onset-aligned spike rates
                alignment = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags) + ...
                    isi;
                alignment_flags = ...
                    valid_time >= alignment + roi_win.(epoch)(1) & ...
                    valid_time < alignment + t2(trial_flags);
                chunk_flags = ...
                    valid_time >= alignment + roi_win.(epoch)(1) & ...
                    valid_time < alignment + roi_win.(epoch)(2);
                aligned_spkrates = spike_rates;
                aligned_spkrates(~alignment_flags') = nan;
                aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                    [roi_length.(epoch),n_flagged_trials])';
                
                % handle cross-validation
                test_idcs = find(~ismember(...
                    flagged_trials,vertcat(xval_train_trials{:,kk})));
                if isempty(test_idcs)
                    continue;
                end
                
                % store tensor & concatenation data
                rand_idcs = randsample(test_idcs,n_concatspercond,true);
                r = aligned_spkrates(rand_idcs,:);
                R(nn,kk+conditions.train.n) = nanmean(r,[1,2]);
            end
        end
        
        % compute scale factors
        scale_factors(:,:,rr,ee) = ...
            R(:,(1:conditions.test.n)+conditions.train.n) ./ R(:,1);
    end
    
    %% plot response dilation
    
    % set dilation boundaries
    dilation_bounds = [-1,1] * 30;
    
    % figure initialization
    fig = figure(figopt,...
        ...'position',[250,730,420,315],...
        'name',sprintf('avgfr_dilation_%s_%s',epoch,contrast_str));
    axes(axesopt.default,...
        'xlim',[1,n_contrasts] + [-1,1]*.65,...
        'xtick',1:n_contrasts,...
        'xticklabel',num2cell(contrast_set),...
        'ylim',dilation_bounds+[-1,1]*range(dilation_bounds)*.05,...
        'ytick',linspace(dilation_bounds(1),dilation_bounds(2),5),...
        'xticklabelrotation',45,...
        ...'ticklength',axesopt.default.ticklength/.75,...
        'layer','bottom',...
        'plotboxaspectratio',[1,1.5,1]);
    title(capitalize(epoch));
    xlabel(sprintf('%s (%s)',contrast_lbl,contrast_units));
    ylabel('<Firing rate> dilation (%)');
    
    % convert scaling factors to dilations
    all_dilations = (nanmean(scale_factors(:,:,:,ee),3) - 1) * 100;
%     all_dilations = nanmean(scale_factors(:,:,:,ee),3) * 3e3;
%     all_dilations = R(:,2:end) * 1.5e0;
%     all_dilations = squeeze(nanmean(zpsths(time>=-334&time<=0,:,:),1)) * 30;
%     all_dilations = squeeze(nanmean(zpsths(time>=0&time<=334,:,:),1)) * 30;
    ref_dilations = all_dilations(:,contrast_mode_idx);
    
    % zero line
    plot(xlim,[1,1]*0,...
        'color','k',...
        'linestyle',':');
    
    % graphical object preallocation
    p = gobjects(3,n_contrasts);
    
    % bar settings
    barwidth = 2/3;
    
    % dilation selection settings
    lowerbound = min(ylim);
    upperbound = max(ylim);
    % lowerbound = -inf;
    % upperbound = +inf;
    
    % average function selection
    avgfun = @(x,d) nanmedian(x,d);
    errfun = @(x) diff(quantile(x,3));
%     avgfun = @(x,d) nanmean(x,d);
%     errfun = @(x) [1,1] * nanstd(x) / sqrt(numel(x));
    
    % compute the reference average
    fov_flags = all(...
        all_dilations >= lowerbound & ...
        all_dilations <= upperbound,2);
    offset = 0; %avgfun(ref_dilations,1); % avgfun(dilations(fov_flags,:),[1,2]);

    %
    dilation_flags = ...
        ...fov_flags & ...
        true(n_neurons,1);
    n_flaggeddilations = sum(dilation_flags)
    dilation_edges = linspace(...
        dilation_bounds(1)*1.25,dilation_bounds(2)*1.25,n_flaggeddilations);
    
    % iterate through conditions
    for kk = 1 : n_contrasts
        
        % mean & standard error across animals
        cond_avg = avgfun(all_dilations(dilation_flags,kk),1) - offset;
        cond_err = errfun(all_dilations(dilation_flags,kk));
        
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
            all_dilations(dilation_flags,kk) - offset,dilation_edges);
        dilation_pdf = 4/15 * barwidth * ...
            (dilation_pdf - min(dilation_pdf)) ./ range(dilation_pdf);
        noise = randn(n_flaggeddilations,1) .* dilation_pdf';
        
        % plot animal threshold dilation
        grapeplot(kk+noise,...
            sort(all_dilations(dilation_flags,kk)) - offset,...
            'marker','o',...
            'markersize',4,...
            'markeredgecolor',contrast_clrs(kk,:),...
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
    x = contrast_set;
    X = [ones(n_contrasts,1),x];
    
    % iterate through units
    for uu = 1 : n_neurons
        
        % linear regression
        thetas = X \ (all_dilations(uu,:)' - offset);
        stretch(uu) = thetas(end);
    end
    
    % one-sample t-test
    [~,pvalue,~,stats.ttest] = ttest(stretch);
    stats.ttest.pvalue = pvalue;
    
    %% plot response stretch
    stretch_bounds = [-1,1] * 1 / (10^strcmpi(contrast_str,'t1'));
    count_bounds = [0,75];
    
    % figure initialization
    fig = figure(figopt,...
        'name',sprintf('avgfr_stretch_%s_%s',epoch,contrast_str));
    axes(...
        axesopt.default,...
        'plotboxaspectratio',[1,3,1],...
        'xlim',count_bounds,...
        'xtick',count_bounds,...
        'xdir','reverse',...
        'ylim',stretch_bounds+[-1,1]*.05*range(stretch_bounds),...
        'ytick',linspace(stretch_bounds(1),stretch_bounds(2),5),...
        'yaxislocation','right');
    title(capitalize(epoch));
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
    text(max(xlim)*.9,.025*range(ylim)+stretch_avg,str,...
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
end