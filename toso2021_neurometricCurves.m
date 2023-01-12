%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% training & test set conditions

% training set conditions
conditions.train = intersectconditions(...
    't1',t1(valid_flags),t_set,...([1:t2_mode_idx-1,t2_mode_idx+1:n_t]),...t_set(t1_mode_idx),...
    'i1',i1(valid_flags),[],...i_set(i1_mode_idx),...
    't2',t2(valid_flags),t_set,...
    'i2',i2(valid_flags),i_set(i2_mode_idx));

% test set conditions
conditions.test = intersectconditions(...
    't1',t1(valid_flags),[],...
    'i1',i1(valid_flags),[],...
    't2',t2(valid_flags),t_set,...
    'i2',i2(valid_flags),i_set);

conditions.train.values
conditions.test.values

%% run settings
n_runs = 5;

%% concatenation settings
n_concatspercond = 2^7; % 2^8
n_concats = n_concatspercond * (conditions.train.n + conditions.test.n);

%% neurometric curve settings
spkintegration_window = min(t_set);

% preallocation
neurocurves = struct();

% iterate through runs
for rr = 1 : n_runs

    %% construct spike counts tensor (neurons X concatenations)

    % preallocation
    concat_spkrates = nan(n_neurons,n_concats);
    concat_t1 = nan(n_concats,1);
    concat_i1 = nan(n_concats,1);
    concat_t2 = nan(n_concats,1);
    concat_i2 = nan(n_concats,1);
    concat_choices = nan(n_concats,1);
    concat_evalset = categorical(nan(n_concats,1),[0,1],{'train','test'});

    % iterate through units
    for nn = 1 : n_neurons
        progressreport(nn,n_neurons,...
            sprintf('sampling concatenations (run %i/%i)',rr,n_runs));
        neuron_flags = data.NeuronNumb == flagged_neurons(nn);

        % preallocation
        xval_train_trials = cell(conditions.test.n,1);

        % iterate through training conditions
        for kk = 1 : conditions.train.n

            % flag trials for the current condition
            t1_flags = ismember(t1,conditions.train.values.t1(kk,:));
            i1_flags = ismember(i1,conditions.train.values.i1(kk,:));
            t2_flags = ismember(t2,conditions.train.values.t2(kk,:));
            i2_flags = ismember(i2,conditions.train.values.i2(kk,:));

            % trial selection
            trial_flags = ...
                valid_flags & ...
                neuron_flags & ...
                t1_flags & ...
                i1_flags & ...
                t2_flags & ...
                i2_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end

            % fetch spike counts & compute spike rates
            spike_counts = data.FR(trial_flags,:);
            spike_rates = ...
                conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;

            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi + ...
                t2(trial_flags);
            alignment_flags = ...
                valid_time >= alignment - spkintegration_window & ...
                valid_time < alignment;
            chunk_flags = alignment_flags;
            aligned_spkrates = spike_rates;
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [spkintegration_window,n_flagged_trials])';
            
            % detect any test conditions that overlap with the current training condition
            t1_xval_condition_flags = any(ismember(...
                conditions.test.values.t1,conditions.train.values.t1(kk,:)),2);
            i1_xval_condition_flags = any(ismember(...
                conditions.test.values.i1,conditions.train.values.i1(kk,:)),2);
            t2_xval_condition_flags = any(ismember(...
                conditions.test.values.t2,conditions.train.values.t2(kk,:)),2);
            i2_xval_condition_flags = any(ismember(...
                conditions.test.values.i2,conditions.train.values.i2(kk,:)),2);
            xval_condition_flags = ...
                t1_xval_condition_flags & ...
                i1_xval_condition_flags & ...
                t2_xval_condition_flags & ...
                i2_xval_condition_flags;
            
            % handle cross-validation for all detected test conditions
            if any(xval_condition_flags)
                xval_condition_idcs = find(xval_condition_flags)';

                % iterate through detected test conditions
                for ii = xval_condition_idcs
                    
                    % detect overlap between training & test trials
                    t1_xval_trial_flags = any(ismember(...
                        t1(flagged_trials),conditions.test.values.t1(ii,:)),2);
                    i1_xval_trial_flags = any(ismember(...
                        i1(flagged_trials),conditions.test.values.i1(ii,:)),2);
                    t2_xval_trial_flags = any(ismember(...
                        t2(flagged_trials),conditions.test.values.t2(ii,:)),2);
                    i2_xval_trial_flags = any(ismember(...
                        i2(flagged_trials),conditions.test.values.i2(ii,:)),2);
                    xval_trial_flags = ...
                        t1_xval_trial_flags & ...
                        i1_xval_trial_flags & ...
                        t2_xval_trial_flags & ...
                        i2_xval_trial_flags;
                    
                    % split the conflicting trials into training & test subsets
                    xval_trials = flagged_trials(xval_trial_flags);
                    n_xval_trials = numel(xval_trials);
                    n_train_trials = round(n_xval_trials * 1 / 2);
                    xval_train_idcs = randperm(n_xval_trials,n_train_trials);
                    xval_train_trials{ii} = xval_trials(xval_train_idcs);
                end
                
                % concatenate sub-sampled training sets across test conditions
                train_idcs = find(ismember(...
                    flagged_trials,vertcat(xval_train_trials{:})));
            else
                
                % train using all trials in the remaining conditions
                train_idcs = 1 : n_flagged_trials;
            end

            % store tensor & concatenation data
            rand_idcs = randsample(train_idcs,n_concatspercond,true);
            concat_idcs = (1 : n_concatspercond) + ...
                n_concatspercond * (kk - 1);
            concat_spkrates(nn,concat_idcs) = ...
                nanmean(aligned_spkrates(rand_idcs,:),2);
            concat_t1(concat_idcs) = t1(flagged_trials(rand_idcs));
            concat_i1(concat_idcs) = i1(flagged_trials(rand_idcs));
            concat_t2(concat_idcs) = t2(flagged_trials(rand_idcs));
            concat_i2(concat_idcs) = i2(flagged_trials(rand_idcs));
            concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
            concat_evalset(concat_idcs) = 'train';
        end

        % iterate through conditions
        for kk = 1 : conditions.test.n

            % flag trials for the current condition
            t1_flags = ismember(t1,conditions.test.values.t1(kk,:));
            i1_flags = ismember(i1,conditions.test.values.i1(kk,:));
            t2_flags = ismember(t2,conditions.test.values.t2(kk,:));
            i2_flags = ismember(i2,conditions.test.values.i2(kk,:));

            % trial selection
            trial_flags = ...
                valid_flags & ...
                neuron_flags & ...
                t1_flags & ...
                i1_flags & ...
                t2_flags & ...
                i2_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end

            % fetch spike counts & compute spike rates
            spike_counts = data.FR(trial_flags,:);
            spike_rates = ...
                conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;

            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi + ...
                t2(trial_flags);
            alignment_flags = ...
                valid_time >= alignment - spkintegration_window & ...
                valid_time < alignment;
            chunk_flags = alignment_flags;
            aligned_spkrates = spike_rates;
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [spkintegration_window,n_flagged_trials])';

            % handle cross-validation
            xval_test_flags = ~ismember(flagged_trials,xval_train_trials{kk});
            test_idcs = find(xval_test_flags);
            if isempty(test_idcs)
                continue;
            end

            % store tensor & concatenation data
            rand_idcs = randsample(test_idcs,n_concatspercond,true);
            concat_idcs = (1 : n_concatspercond) + ...
                n_concatspercond * (kk + conditions.train.n - 1);
            concat_spkrates(nn,concat_idcs) = ...
                nanmean(aligned_spkrates(rand_idcs,:),2);
            concat_t1(concat_idcs) = t1(flagged_trials(rand_idcs));
            concat_i1(concat_idcs) = i1(flagged_trials(rand_idcs));
            concat_t2(concat_idcs) = t2(flagged_trials(rand_idcs));
            concat_i2(concat_idcs) = i2(flagged_trials(rand_idcs));
            concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
            concat_evalset(concat_idcs) = 'test';
        end
    end

    %% construct neurometric curves

    % flag training concatenations
    train_flags = concat_evalset == 'train';

    % linear discriminant analysis
    X = concat_spkrates(:,train_flags)';
    threshold = median(s1(valid_flags));
    y = concat_t2(train_flags) > threshold;
    within_class_var = cell2mat(...
        arrayfun(@(x)nanvar(X(y==x,:)),unique(y),...
        'uniformoutput',false));
    invalid_flags = any(within_class_var == 0 | isnan(within_class_var));
    lda_mdl = fitcdiscr(X(:,~invalid_flags),y,...
        'discrimtype','linear');

    % flag test concatenations
    test_flags = concat_evalset == 'test';

    % preallocation
    neuro_choices = nan(n_concats,1);

    % neural "judgments"
    neuro_choices(test_flags) = ...
        lda_mdl.predict(concat_spkrates(~invalid_flags,test_flags)');

    %% construct neurophysical triple

    % iterate through contrasts
    for kk = 1 : n_contrasts
        contrast_flags = concat_i2 == contrast_set(kk);

        % iterate through stimuli
        for ii = 1 : n_stimuli
            stimulus_flags = concat_t2 == stim_set(ii);
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

%% visualize population state at T2 offset
X = concat_spkrates(~invalid_flags,:)';

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
    stimulus_flags = concat_t2 == stim_set(ii);

    % plot state space projections
    plot(score(stimulus_flags,1),...
        score(stimulus_flags,2),...
        'linestyle','none',...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6,...
        'markerfacecolor','none',...
        'markeredgecolor','k');
end

% iterate through stimuli
for ii = 1 : n_stimuli
    stimulus_flags = concat_t2 == stim_set(ii);

    % plot state space projections
    plot(score(stimulus_flags,1),...
        score(stimulus_flags,2),...
        'linestyle','none',...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6-1.5,...
        'markerfacecolor',s2_clrs(ii,:),...
        'markeredgecolor','none');
%     plot3(score(stimulus_flags,1),...
%         score(stimulus_flags,2),...
%         score(stimulus_flags,3),...
%         'linestyle','none',...
%         'linewidth',1.5,...
%         'marker','o',...
%         'markersize',6-1.5,...
%         'markerfacecolor',s2_clrs(ii,:),...
%         'markeredgecolor','none');
end

% update axes
xlim(xlim + [-1,1] * .05 * range(xlim));
ylim(ylim + [-1,1] * .05 * range(ylim));

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

% iterate through stimuli
for ii = 1 : n_stimuli
    stimulus_flags = concat_t2 == stim_set(ii);

    % plot projections onto linear discriminant
    bincounts = histcounts(score(stimulus_flags),binedges);
    h = histogram(...
        'binedges',binedges,...
        'bincounts',bincounts + prevcounts,...
        'linewidth',1.5,...
        'edgecolor','none',...
        'facecolor',s2_clrs(ii,:),...
        'facealpha',1);

    % update running counts
    prevcounts = prevcounts + bincounts;
    if stim_set(ii) < threshold
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
psyopt.fit.borders = [0,1;0,1;0,.25;0,.25;0,0];
psyopt.fit.fixedPars = [nan,nan,nan,nan,0];
psyopt.fit.stepN = [100,100,20,20,20];

% iterate through contrasts
for kk = 1 : n_contrasts

    % pool neural choices across runs
    neurocurve_pools(kk).x = unique(neurocurves(kk).x);
    neurocurve_pools(kk).y = sum(neurocurves(kk).y,2);
    neurocurve_pools(kk).n = sum(neurocurves(kk).n,2);
    neurocurve_pools(kk).err = ...
        std(neurocurves(kk).y ./ neurocurves(kk).n,0,2) / sqrt(n_runs);

    % neurophysical triple
    neuro_triple = [...
        neurocurve_pools(kk).x,...
        neurocurve_pools(kk).y,...
        neurocurve_pools(kk).n...
        ];
%     if kk ~= 2
%         neuro_triple(6,2) = nan;
%     end
% neuro_triple

    % fit psychometric curve
    neurocurve_pools(kk).fit = psignifit(neuro_triple,psyopt.fit);
end

%% plot neurometric function

% figure initialization
fig = figure(figopt,...
    'name',sprintf('neurometric_curves_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
title('Neurometric curves');
xlabel(sprintf('%.2f \\times T_2 + %.2f \\times T_1 (ms)',w2,w1));
ylabel('P(T_2 > T_1)');

% psychometric plot settings
psyopt.plot = struct();
psyopt.plot.linewidth = 1.5;
psyopt.plot.marker = 's';
psyopt.plot.markersize = 9;
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
    psyopt.plot.normalizemarkersize = true;
    psyopt.plot.plotfit = sum(neurocurve_pools(kk).n ~= 0) >= 2;
    p(kk) = plotpsy(neurocurve_pools(kk),neurocurve_pools(kk).fit,psyopt.plot);
end

%% inset with delta P(long)
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.inset.se,...
    'ylim',[-.25,.25]+[-1,1]*.05*.75,...
    'ytick',[-.25,0,.25],...
    'clipping','off');
ylabel('\DeltaP(T_2 > T_1)',...
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
    p_ctrl = neurocurve_pools(contrast_mode_idx).y ./ ...
        neurocurve_pools(contrast_mode_idx).n;
    p_cond = neurocurve_pools(kk).y ./ neurocurve_pools(kk).n;
    delta_p = p_cond - p_ctrl;
    err = sqrt((neurocurve_pools(kk).err .^ 2) + (neurocurve_pools(kk).err .^ 2));
    errorbar(normstim_set,delta_p,err,...
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
