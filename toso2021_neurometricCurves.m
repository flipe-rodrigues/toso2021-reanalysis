%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% contrast settings
contrast_str = 'i2';
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str(1:end-1),'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);

%% stimulus selection (for training & test sets)

% stimuli
stim2train_idcs = [1 : t2_mode_idx - 1, t2_mode_idx + 1 : n_stimuli];
stim2train_n = numel(stim2train_idcs);
stim2test_idcs = 1 : n_stimuli;
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

%% neurometric curve settings
spkintegration_window = 100;
n_runs = 25;

% preallocation
neurocurves = struct();

% stimulus 1 flags
i1_flags = i1 == i_set(i1_mode_idx);
t1_flags = t1 > t_set(2) & t1 < t_set(end-1);

% iterate through runs
for rr = 1 : n_runs
    
    %% construct spike rate tensor (time X neurons X concatenations)
    
    % preallocation
    concat_spkrates = nan(n_neurons,n_concats);
    concat_contrasts = nan(n_concats,1);
    concat_stimuli = nan(n_concats,1);
    concat_choices = nan(n_concats,1);
    concat_evalset = nan(n_concats,1);
    
    % iterate through units
    for nn = 1 : n_neurons
        progressreport(nn,n_neurons,...
            sprintf('constructing concatenations on run %i',rr));
        neuron_flags = data.NeuronNumb == flagged_neurons(nn);
        
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
                i1_flags & ...
                t1_flags & ...
                stimulus_flags & ...
                contrast_flags;
            flagged_trials = find(s2_spike_flags);
            if sum(s2_spike_flags) == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            s2_spike_counts = data.FR(s2_spike_flags,:);
            s2_spike_rates = ...
                conv2(1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
            n_trials = size(s2_spike_counts,1);
            
            % T2-offset-aligned spike rates
            s2_alignment_offset = ...
                pre_init_padding + ...
                pre_t1_delay(s2_spike_flags) + ...
                t1(s2_spike_flags) + ...
                inter_t1t2_delay + ...
                t2(s2_spike_flags);
            s2_alignment_flags = ...
                valid_time >= s2_alignment_offset - spkintegration_window & ...
                valid_time < s2_alignment_offset;
            s2_chunk_flags = s2_alignment_flags;
            s2_spkrates = s2_spike_rates;
            s2_spkrates(~s2_alignment_flags') = nan;
            s2_spkrates = ...
                reshape(s2_spkrates(s2_chunk_flags'),[spkintegration_window,n_trials])';
            
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
            concat_spkrates(nn,concat_idcs) = ...
                nanmean(s2_spkrates(rand_idcs,:),2);
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
                continue;
            end
            
            % fetch spike counts & compute spike rates
            s2_spike_counts = data.FR(s2_spike_flags,:);
            s2_spike_rates = ...
                conv2(1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
            n_trials = size(s2_spike_counts,1);
            
            % T2-offset-aligned spike rates
            s2_alignment_offset = ...
                pre_init_padding + ...
                pre_t1_delay(s2_spike_flags) + ...
                t1(s2_spike_flags) + ...
                inter_t1t2_delay + ...
                t2(s2_spike_flags);
            s2_alignment_flags = ...
                valid_time >= s2_alignment_offset - spkintegration_window & ...
                valid_time < s2_alignment_offset;
            s2_chunk_flags = s2_alignment_flags;
            s2_spkrates = s2_spike_rates;
            s2_spkrates(~s2_alignment_flags') = nan;
            s2_spkrates = ...
                reshape(s2_spkrates(s2_chunk_flags'),[spkintegration_window,n_trials])';
            
            % handle cross-validation
            test_flags = ~ismember(flagged_trials,nn_train_trials);
            test_idcs = find(test_flags);
            
            % store tensor & concatenation data
            rand_idcs = randsample(test_idcs,n_concatspercond,true);
            concat_idcs = (1 : n_concatspercond) + ...
                n_concatspercond * (kk + conditions.train.n - 1);
            concat_spkrates(nn,concat_idcs) = ...
                nanmean(s2_spkrates(rand_idcs,:),2);
            concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
            concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
            concat_choices(concat_idcs) = choices(flagged_trials(rand_idcs));
            concat_evalset(concat_idcs) = 1;
        end
    end
    
    % convert train / test flags to categorical
    concat_evalset = categorical(concat_evalset,[0,1],{'train','test'});
    
    %% construct neurometric curves
    
    % flag training concatenations
    train_flags = concat_evalset == 'train';
    train_idcs = find(train_flags);
    
    % linear discriminant analysis
    X = concat_spkrates(:,train_flags)';
    y = concat_stimuli(train_flags) > median(t1);
    lda_mdl = fitcdiscr(X,y,...
        'discrimtype','linear');
    
    % flag test concatenations
    test_flags = concat_evalset == 'test';
    
    % preallocation
    neuro_choices = nan(n_concats,1);

    % neural "judgements"
    neuro_choices(test_flags) = ...
        lda_mdl.predict(concat_spkrates(:,test_flags)');
    
    %% construct neurophysical triple
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        contrast_flags = concat_contrasts == contrast_set(kk);
        
        % iterate through stimuli
        for ii = 1 : n_stimuli
            stimulus_flags = concat_stimuli == stim_set(ii);
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
X = concat_spkrates';

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
    stimulus_flags = concat_stimuli == stim_set(ii);
    
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
    stimulus_flags = concat_stimuli == stim_set(ii);
    
    % plot state space projections
    plot(score(stimulus_flags,1),...
        score(stimulus_flags,2),...
        'linestyle','none',...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6-1.5,...
        'markerfacecolor',t2_clrs(ii,:),...
        'markeredgecolor','none');
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
shortcounts = zeros(1,n_bins);
longcounts = zeros(1,n_bins);

% iterate through stimuli
for ii = 1 : n_stimuli
    stimulus_flags = concat_stimuli == stim_set(ii);
    
    % plot projections onto linear discriminant
    bincounts = histcounts(score(stimulus_flags),binedges);
    h = histogram(...
        'binedges',binedges,...
        'bincounts',bincounts + prevcounts,...
        'linewidth',1.5,...
        'edgecolor','none',...
        'facecolor',t2_clrs(ii,:),...
        'facealpha',1);
    
    % update running counts
    prevcounts = prevcounts + bincounts;
    if stim_set(ii) < median(t1)
        shortcounts = shortcounts + bincounts;
    end
    if stim_set(ii) > median(t1)
        longcounts = longcounts + bincounts;
    end
    
    % manage ui stack
    uistack(h,'bottom');
end

% plot distribution outline
stairs(binedges,[prevcounts,prevcounts(end)],...
    'linewidth',1.5,...
    'color','k');
stairs(binedges,[shortcounts,shortcounts(end)],...
    'linewidth',1.5,...
    'linestyle',':',...
    'color','k');
stairs(binedges,[longcounts,longcounts(end)],...
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