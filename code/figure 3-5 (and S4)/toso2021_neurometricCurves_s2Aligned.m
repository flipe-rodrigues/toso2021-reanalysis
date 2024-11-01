%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% run settings
n_runs = 100;

%% concatenation settings
n_concats_max = 2^7;

%% training & test set conditions

% training set conditions
if strcmpi(contrast_str,'t1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'t2')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set,[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'choice')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'correct')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
end

% test set conditions
if strcmpi(contrast_str,'t1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'t2')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set,[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'choice')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),choice_set,[],...
        'correct',correct(valid_flags),[],[]);
elseif strcmpi(contrast_str,'correct')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),t_set,[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),correct_set,[]);
end

% print training & test conditions
fprintf('\nTRAINING CONDITIONS:\n');
disp(conditions.train.values);
fprintf('\nTEST CONDITIONS:\n');
conditions.test.values

%% subject selection
subject_flags = ismember(subjects,subject_set);

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
n_concats_train = round(conditions.train.weights * n_concats_max);
n_concats_test = round(conditions.test.weights * n_concats_max);
n_concats_total = sum([n_concats_train; n_concats_test]);

%% neurometric curve settings
spk_integration_win = min(t_set);
n_integration_bins = spk_integration_win / psthbin;

% preallocation
neurocurves = struct();

% selection criteria
n_trials_cutoff = 0;

% iterate through runs
for rr = 1 : n_runs
    
    %% construct spike counts tensor (neurons X concatenations)
    
    % preallocation
    concat_spkcounts = nan(n_neurons,n_concats_total);
    concat_s1 = nan(n_concats_total,1);
    concat_s2 = nan(n_concats_total,1);
    concat_contrasts = nan(n_concats_total,1);
    concat_choices = nan(n_concats_total,1);
    concat_evalset = ...
        categorical(nan(n_concats_total,1),[0,1],{'train','test'});
    
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
                subject_flags & ...
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials <= n_trials_cutoff
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_counts = data.SDF(trial_flags,:)';
            
            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi + ...
                t2(trial_flags);
            alignment_flags = ...
                padded_time >= alignment - spk_integration_win & ...
                padded_time < alignment;
            chunk_flags = alignment_flags;
            aligned_spkcounts = spike_counts;
            aligned_spkcounts(~alignment_flags') = nan;
            aligned_spkcounts = reshape(aligned_spkcounts(chunk_flags'),...
                [n_integration_bins,n_flagged_trials])';
            
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
            rand_idcs = randsample(train_idcs,n_concats_train(kk),true);
            concat_idcs = (1 : n_concats_train(kk)) + ...
                sum(n_concats_train(1:kk-1));
            concat_spkcounts(nn,concat_idcs) = ...
                nanmean(aligned_spkcounts(rand_idcs,:),2);
            concat_s1(concat_idcs) = s1(flagged_trials(rand_idcs));
            concat_s2(concat_idcs) = s2(flagged_trials(rand_idcs));
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
                subject_flags & ...
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials <= n_trials_cutoff
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_counts = data.SDF(trial_flags,:)';
            
            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi + ...
                t2(trial_flags);
            alignment_flags = ...
                padded_time >= alignment - spk_integration_win & ...
                padded_time < alignment;
            chunk_flags = alignment_flags;
            aligned_spkcounts = spike_counts;
            aligned_spkcounts(~alignment_flags') = nan;
            aligned_spkcounts = reshape(aligned_spkcounts(chunk_flags'),...
                [n_integration_bins,n_flagged_trials])';
            
            % handle cross-validation
            test_idcs = find(~ismember(...
                flagged_trials,vertcat(xval_train_trials{:,kk})));
            if isempty(test_idcs)
                continue;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(test_idcs,n_concats_test(kk),true);
            concat_idcs = (1 : n_concats_test(kk)) + ...
                sum(n_concats_test(1:kk-1)) + ...
                sum(n_concats_train);
            concat_spkcounts(nn,concat_idcs) = ...
                nanmean(aligned_spkcounts(rand_idcs,:),2);
            concat_s1(concat_idcs) = s1(flagged_trials(rand_idcs));
            concat_s2(concat_idcs) = s2(flagged_trials(rand_idcs));
            concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
            concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
            concat_evalset(concat_idcs) = 'test';
        end
    end
    
    %% construct neurometric curves
    
    % flag training concatenations
    train_flags = concat_evalset == 'train';
    
    % linear discriminant analysis
    X = concat_spkcounts(:,train_flags)';
    y = concat_s2(train_flags) > concat_s1(train_flags);
    within_class_var = cell2mat(...
        arrayfun(@(x)nanvar(X(y==x,:)),unique(y),...
        'uniformoutput',false));
    invalid_flags = ...
        any(within_class_var == 0 | isnan(within_class_var));
    lda_mdl = fitcdiscr(X(:,~invalid_flags),y,...
        'discrimtype','linear');
    
    % flag test concatenations
    test_flags = concat_evalset == 'test';
    
    % preallocation
    neuro_choices = nan(n_concats_total,1);
    
    % neural "judgments"
    neuro_choices(test_flags) = ...
        lda_mdl.predict(concat_spkcounts(~invalid_flags,test_flags)');
    
    %% construct neurophysical triple
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        contrast_flags = concat_contrasts == contrast_set(kk);
        
        % iterate through stimuli
        for ii = 1 : n_stimuli
            stimulus_flags = concat_s2 == stim_set(ii);
            concat_flags = ...
                test_flags & ...
                contrast_flags & ...
                stimulus_flags;
            neurocurves(kk).x(ii,rr) = normstim_set(ii);
            neurocurves(kk).y(ii,rr) = sum(neuro_choices(concat_flags),[1,2]);
            neurocurves(kk).n(ii,rr) = sum(concat_flags);
            neurocurves(kk).err(ii,rr) = ...
                std(neuro_choices(concat_flags),0,[1,2]) / sqrt(sum(concat_flags));
        end
    end
end

%% visualize population state at S2 offset
X = concat_spkcounts(:,train_flags)';

% normalization
Z = zscore(X);

% pca
[coeff,score,~,~,explained] = pca(X);

% figure initialization
fig = figure(figopt,...
    'name','pca_visualization');

% axes initialization
set(gca,...
    axesopt.default,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'xtick',0,...
    'ytick',0);
xlabel(sprintf('%s\n%.1f%% variance','PC 1',explained(1)),...
    'horizontalalignment','center');
ylabel(sprintf('%s\n%.1f%% variance','PC 2',explained(2)),...
    'horizontalalignment','center');

% iterate through stimuli
for ii = 1 : n_stimuli
    stimulus_flags = concat_s2 == stim_set(ii);
    concat_flags = ...
        train_flags & ...
        stimulus_flags;
    
    % plot state space projections
    plot(score(concat_flags,1),...
        score(concat_flags,2),...
        'linestyle','none',...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6,...
        'markerfacecolor','none',...
        'markeredgecolor','k');
end

% iterate through stimuli
for ii = 1 : n_stimuli
    stimulus_flags = concat_s2 == stim_set(ii);
    concat_flags = ...
        train_flags & ...
        stimulus_flags;
    
    % plot state space projections
    plot(score(concat_flags,1),...
        score(concat_flags,2),...
        'linestyle','none',...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6-1.5,...
        'markerfacecolor',s2_clrs(ii,:),...
        'markeredgecolor','none');
end

% update axes
xlim(xlim + [-1,1] * .05 * range(xlim));
ylim(ylim + [-1,1] * .05 * range(ylim));

% plot decision boundary
betas = lda_mdl.Coeffs(2,1).Linear' * coeff(:,1:2);
x = linspace(min(xlim),max(xlim),1e3);
y = betas(1) + betas(2) * x;
plot(x,y,'--k',...
    'linewidth',1.5);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% LDA visualization

% figure initialization
fig = figure(figopt,...
    'name','lda_visualization');

% axes initialization
set(gca,...
    axesopt.default,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'xtick',[],...
    'ytick',[],...
    'ycolor','none',...
    'plotboxaspectratio',[3,1,1]);
xlabel('Projection onto linear discriminant');

% project data onto linear discriminant
score = X * lda_mdl.Coeffs(2,1).Linear;

% bin settings
n_bins = 50;
binedges = linspace(min(score)*.95,max(score)*1.05,n_bins+1);
prevcounts = zeros(1,n_bins);
s2lessthans1_counts = zeros(1,n_bins);
s2greaterthans1_counts = zeros(1,n_bins);

%
threshold = nanmedian(stimuli(valid_flags));

% iterate through stimuli
for ii = 1 : n_stimuli
    concat_flags = concat_s2(train_flags) == stim_set(ii);
    
    % plot projections onto linear discriminant
    bincounts = histcounts(score(concat_flags),binedges);
    h = histogram(...
        'binedges',binedges,...
        'bincounts',bincounts + prevcounts,...
        'linewidth',1.5,...
        'edgecolor','none',...
        'facecolor',s2_clrs(ii,:),...
        'facealpha',1);
    
    % update running counts
    prevcounts = prevcounts + bincounts;
    if stim_set(ii) <= threshold
        s2lessthans1_counts = s2lessthans1_counts + bincounts;
    end
    if stim_set(ii) > threshold
        s2greaterthans1_counts = s2greaterthans1_counts + bincounts;
    end
    
    % manage ui stack
    uistack(h,'bottom');
end

% plot distribution outline
stairs(binedges,[prevcounts,prevcounts(end)],...
    'linewidth',1.5,...
    'color','k');
stairs(binedges,[s2lessthans1_counts,s2lessthans1_counts(end)],...
    'linewidth',1.5,...
    'linestyle',':',...
    'color','k');
stairs(binedges,[s2greaterthans1_counts,s2greaterthans1_counts(end)],...
    'linewidth',1.5,...
    'linestyle',':',...
    'color','k');

% plot decision boundary
plot([1,1]*mean(xlim),ylim,'--k',...
    'linewidth',1.5);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% fit neurometric function

% preallocation
neurocurve_pools = struct();

% psychometric fit settings
psyopt.fit = struct();
psyopt.fit.expType = 'YesNo';
psyopt.fit.sigmoidName = 'rgumbel';
psyopt.fit.estimateType = 'MAP';
psyopt.fit.confP = [.95,.9,.68];
psyopt.fit.borders = [0,1; 0,1; 0,.25; 0,.25; 0,0];
psyopt.fit.fixedPars = [nan,nan,nan,nan,0];
psyopt.fit.stepN = [100,100,40,40,20];

% average & error function selection
avgfun = @(x,d) nansum(x,d);
errfun = @(x,d) std(x,0,d) / 1 .* [1,1];

% iterate through contrasts
for kk = 1 : n_contrasts
    
    % pool neural choices across runs
    neurocurve_pools(kk).x = unique(neurocurves(kk).x);
    neurocurve_pools(kk).y = avgfun(neurocurves(kk).y,2);
    neurocurve_pools(kk).n = avgfun(neurocurves(kk).n,2);
    neurocurve_pools(kk).err = ....
        errfun(neurocurves(kk).y ./ neurocurves(kk).n,2);
    
    % neurophysical triple
    neuro_triple = [...
        neurocurve_pools(kk).x,...
        neurocurve_pools(kk).y,...
        neurocurve_pools(kk).n...
        ];
    
    % fit psychometric curve
    neurocurve_pools(kk).fit = psignifit(neuro_triple,psyopt.fit);
end

%% plot neurometric function

% figure initialization
fig = figure(figopt,...
    'name',sprintf('neurometric_curves_s2_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
xlabel(sprintf('%s (%s)',s2_lbl,s_units));
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% psychometric plot settings
psyopt.plot = struct();
psyopt.plot.linewidth = 1.5;
psyopt.plot.marker = 's';
psyopt.plot.markersize = 9.5;
psyopt.plot.plotdata = true;
psyopt.plot.gradeclrs = false;
psyopt.plot.patchci = false;
psyopt.plot.overallvisibility = 'off';
psyopt.plot.normalizemarkersize = true;
psyopt.plot.plotfit = true;

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*median(normstim_set),ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through contrasts
for kk = 1 : n_contrasts
    
    % plot psychometric curve
    psyopt.plot.datafaceclr = contrast_clrs(kk,:);
    psyopt.plot.overallvisibility = 'off';
    psyopt.plot.plotfit = sum(neurocurve_pools(kk).n ~= 0) >= 2;
    p(kk) = plotpsy(neurocurve_pools(kk),neurocurve_pools(kk).fit,psyopt.plot);
end

%% inset with delta P(S2 > S1)
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.inset.se,...
    'ylim',[-.25,.25]+[-1,1]*.05*.75,...
    'ytick',[-.25,0,.25],...
    'clipping','off');
ylabel(sprintf('\DeltaP(%s > %s)',s2_lbl,s1_lbl),...
    'rotation',-90,...
    'verticalalignment','bottom');

% categorical boundary
plot([1,1]*median(normstim_set),ylim,...
    'color','k',...
    'linestyle',':');

% zero level
plot(xlim,[0,0],...
    'color','k',...
    'linestyle',':');

% iterate through contrasts
for kk = 1 : n_contrasts
    if kk == contrast_mode_idx
        continue;
    end
    
    % plot delta P(S2 > S1)
    p_ctrl = neurocurve_pools(contrast_mode_idx).y ./ ...
        neurocurve_pools(contrast_mode_idx).n;
    p_cond = neurocurve_pools(kk).y ./ neurocurve_pools(kk).n;
    delta_p = p_cond - p_ctrl;
    err = sqrt((neurocurve_pools(kk).err .^ 2) + (neurocurve_pools(kk).err .^ 2));
    errorbar(normstim_set,delta_p,err(:,1),err(:,2),...
        'color','k',...
        'marker',psyopt.plot.marker,...
        'markersize',psyopt.plot.markersize-2.5,...
        'markerfacecolor',contrast_clrs(kk,:),...
        'linewidth',psyopt.plot.linewidth,...
        'capsize',0);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end