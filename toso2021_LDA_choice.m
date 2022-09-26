%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% concatenation settings
n_concatsperchoice = 2^10;
n_concats = n_concatsperchoice * n_choices;

%% LDA settings
spkintegration_window = 100;

%% construct concatenations for LDA

% preallocation
concat_spkrates = nan(n_neurons,n_concats);
concat_choices = nan(n_concats,1);

% iterate through units
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing concatenations for LDA');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % preallocation
    nn_train_trials = cell(conditions.train.n,1);
    
    % iterate through choices
    for cc = 1 : n_choices
        choice_flags = choice == choice_set(cc);
        
        % flag trials for the current condition
        s2_spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            choice_flags;
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
        s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(s2_spike_flags) + ...
            t1(s2_spike_flags) + ...
            isi + ...
            t2(s2_spike_flags);
        s2_alignment_flags = ...
            valid_time >= s2_alignment + post_s2_delay - spkintegration_window & ...
            valid_time < s2_alignment + post_s2_delay;
        s2_chunk_flags = s2_alignment_flags;
        s2_spkrates = s2_spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = ...
            reshape(s2_spkrates(s2_chunk_flags'),[spkintegration_window,n_trials])';
        
        % store tensor & concatenation data
        rand_idcs = randsample(n_trials,n_concatsperchoice,true);
        concat_idcs = (1 : n_concatsperchoice) + n_concatsperchoice * (cc - 1);
        concat_spkrates(nn,concat_idcs) = ...
            nanmean(s2_spkrates(rand_idcs,:),2);
        concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
    end
end

%% LDA

% linear discriminant analysis
X = concat_spkrates';
y = concat_choices;
lda_mdl = fitcdiscr(X,y,...
    'discrimtype','linear');

% fisher's linear discriminant
sig = cov(x0);
mu0 = mean(x0);
mu1 = mean(x1);
w = sig \ (mu0 - mu1)';

%% visualize population state at ROI
X = concat_spkrates';

% normalization
Z = zscore(X);

% pca
[coeff,score,~,~,explained] = pca(X);

% figure initialization
fig = figure(figopt,...
    'name','pca_visualization_choice');

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

% iterate through choices
for cc = 1 : n_choices
    choice_flags = concat_choices == choice_set(cc);
    
    % plot state space projections
    plot(score(choice_flags,1),...
        score(choice_flags,2),...
        ...score(choice_flags,3),...
        'linestyle','none',...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6,...
        'markerfacecolor','none',...
        'markeredgecolor','k');
end

% iterate through choices
for cc = 1 : n_choices
    choice_flags = concat_choices == choice_set(cc);
    
    % plot state space projections
    plot(score(choice_flags,1),...
        score(choice_flags,2),...
        ...score(choice_flags,3),...
        'linestyle','none',...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6-1.5,...
        'markerfacecolor',choice_clrs(cc,:),...
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
    'name','lda_visualization_choice');

% axes initialization
set(gca,...
    axesopt.default,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'xtick',[],...
    'ytick',[],...
    'ycolor','none',...
    'plotboxaspectratio',[3,1,1]);
xlabel('Projection onto linear discriminant (LD)');

% project data onto linear discriminant
score = X * lda_mdl.Coeffs(2,1).Linear;

% bin settings
n_bins = 50;
binedges = linspace(min(score)*.95,max(score)*1.05,n_bins+1);

% iterate through choices
for cc = 1 : n_choices
    choice_flags = concat_choices == choice_set(cc);
    
    % plot projections onto linear discriminant
    bincounts = histcounts(score(choice_flags),binedges);
    h = histogram(...
        'binedges',binedges,...
        'bincounts',bincounts,...
        'linewidth',1.5,...
        'edgecolor','none',...
        'facecolor',choice_clrs(cc,:),...
        'facealpha',1);

    % plot distribution outline
    stairs(binedges,[bincounts,bincounts(end)],...
        'linewidth',1.5,...
        'color','k');
    
    % manage ui stack
    uistack(h,'bottom');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% epoch settings

% preallocation
epochs = struct();
psths = struct();
zpsths = struct();

% epoch alignments
epochs.label.pre_s1 = 'pre-T_{1}-delay';
epochs.label.s1 = 'T_{1}';
epochs.label.isi = 'ISI';
epochs.label.s2 = 'T_{2}';
epochs.label.post_s2 = 'post-T_{2}-delay';
epochs.label.go = 'go cue';

% epoch rois
epochs.roi.pre_s1 = [-500,unique(pre_t1_delay(valid_flags))];
epochs.roi.s1 = [0,t_set(end)];
epochs.roi.isi = [0,isi];
epochs.roi.s2 = [0,t_set(end)];
epochs.roi.post_s2 = [0,1] * post_t2_delay;
epochs.roi.go = [0,450];

% epoch contrasts
epochs.contrast.pre_s1 = 't1';
epochs.contrast.s1 = 't1';
epochs.contrast.isi = 't1';
epochs.contrast.s2 = 't1';
epochs.contrast.post_s2 = 't1';
epochs.contrast.go = 't1';

% iterate through epochs
epoch_labels = fieldnames(epochs.label);
n_epochs = numel(epoch_labels);
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    epoch_contrast_str = epochs.contrast.(epoch);
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    epochs.n_bins.(epoch) = range(epochs.roi.(epoch)) * psthbin;
    epochs.time.(epoch) = linspace(epochs.roi.(epoch)(1),epochs.roi.(epoch)(2),epochs.n_bins.(epoch));
    psths.(epoch) = nan(epochs.n_bins.(epoch),n_neurons,epoch_n_contrasts);
end

%% for cross-epoch intersectional selection
surviving_neurons = true(n_neurons,1);

%% construct pre-S1-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing pre-s1 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.go;
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % pre T1 delay-aligned spike rates
        pre_s1_alignment = ...
            repmat(pre_init_padding,n_trials,1);
        pre_s1_alignment_flags = ...
            valid_time >= pre_s1_alignment + epochs.roi.pre_s1(1) & ...
            valid_time < pre_s1_alignment + pre_t1_delay(spike_flags);
        pre_s1_chunk_flags = ...
            valid_time >= pre_s1_alignment + epochs.roi.pre_s1(1) & ...
            valid_time < pre_s1_alignment + epochs.roi.pre_s1(2);
        pre_s1_spkrates = spike_rates;
        pre_s1_spkrates(~pre_s1_alignment_flags') = nan;
        pre_s1_spkrates = reshape(...
            pre_s1_spkrates(pre_s1_chunk_flags'),[epochs.n_bins.pre_s1,n_trials])';
        
        % compute mean spike density function
        psths.pre_s1(:,nn,ii) = nanmean(pre_s1_spkrates,1);
    end
end

%% construct S1-onset-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing s1 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.go;
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % T1-aligned spike rates
        s1_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags);
        s1_alignment_flags = ...
            valid_time >= s1_alignment + epochs.roi.s1(1) & ...
            valid_time < s1_alignment + t1(spike_flags);
        s1_chunk_flags = ...
            valid_time >= s1_alignment + epochs.roi.s1(1) & ...
            valid_time < s1_alignment + epochs.roi.s1(2);
        s1_spkrates = spike_rates;
        s1_spkrates(~s1_alignment_flags') = nan;
        s1_spkrates = reshape(...
            s1_spkrates(s1_chunk_flags'),[epochs.n_bins.s1,n_trials])';
        
        % compute mean spike density function
        psths.s1(:,nn,ii) = nanmean(s1_spkrates,1);
    end
end

%% construct ISI-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing isi psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.go;
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % ISI-aligned spike rates
        isi_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags);
        isi_alignment_flags = ...
            valid_time >= isi_alignment + epochs.roi.isi(1) & ...
            valid_time < isi_alignment + isi;
        isi_chunk_flags = ...
            valid_time >= isi_alignment + epochs.roi.isi(1) & ...
            valid_time < isi_alignment + epochs.roi.isi(2);
        isi_spkrates = spike_rates;
        isi_spkrates(~isi_alignment_flags') = nan;
        isi_spkrates = reshape(...
            isi_spkrates(isi_chunk_flags'),[epochs.n_bins.isi,n_trials])';
        
        % compute mean spike density function
        psths.isi(:,nn,ii) = nanmean(isi_spkrates,1);
    end
end

%% construct S2-onset-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing s2 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.go;
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % T2-aligned spike rates
        s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi;
        s2_alignment_flags = ...
            valid_time >= s2_alignment + epochs.roi.s2(1) & ...
            valid_time < s2_alignment + t2(spike_flags);
        s2_chunk_flags = ...
            valid_time >= s2_alignment + epochs.roi.s2(1) & ...
            valid_time < s2_alignment + epochs.roi.s2(2);
        s2_spkrates = spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(...
            s2_spkrates(s2_chunk_flags'),[epochs.n_bins.s2,n_trials])';
        
        % compute mean spike density function
        psths.s2(:,nn,ii) = nanmean(s2_spkrates,1);
    end
end

%% construct post-S2-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing post-s2 psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.go;
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % post T2 delay-aligned spike rates
        post_s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi + ...
            t2(spike_flags);
        post_s2_alignment_flags = ...
            valid_time >= post_s2_alignment + epochs.roi.post_s2(1) & ...
            valid_time < post_s2_alignment + post_t2_delay;
        post_s2_chunk_flags = ...
            valid_time >= post_s2_alignment + epochs.roi.post_s2(1) & ...
            valid_time < post_s2_alignment + epochs.roi.post_s2(2);
        post_s2_spkrates = spike_rates;
        post_s2_spkrates(~post_s2_alignment_flags') = nan;
        post_s2_spkrates = reshape(...
            post_s2_spkrates(post_s2_chunk_flags'),[epochs.n_bins.post_s2,n_trials])';
        
        % compute mean spike density function
        psths.post_s2(:,nn,ii) = nanmean(post_s2_spkrates,1);
    end
end

%% construct go-cue-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'constructing go psths');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % contrast settings
    epoch_contrast_str = epochs.contrast.go;
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    
    % iterate through contrasts
    for ii = 1 : epoch_n_contrasts
        contrast_flags = epoch_contrasts == epoch_contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            surviving_neurons(nn) = false;
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % go-aligned spike rates
        go_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi + ...
            t2(spike_flags) + ...
            post_t2_delay;
        go_alignment_flags = ...
            valid_time >= go_alignment + epochs.roi.go(1) & ...
            valid_time < go_alignment + epochs.roi.go(2);
        go_chunk_flags = go_alignment_flags;
        go_spkrates = spike_rates;
        go_spkrates(~go_alignment_flags') = nan;
        go_spkrates = reshape(...
            go_spkrates(go_chunk_flags'),[epochs.n_bins.go,n_trials])';
        
        % compute mean spike density function
        psths.go(:,nn,ii) = nanmean(go_spkrates,1);
    end
end

%% intersectional neuron selection
psths.pre_s1 = psths.pre_s1(:,surviving_neurons,:);
psths.s1 = psths.s1(:,surviving_neurons,:);
psths.isi = psths.isi(:,surviving_neurons,:);
psths.s2 = psths.s2(:,surviving_neurons,:);
psths.post_s2 = psths.post_s2(:,surviving_neurons,:);
psths.go = psths.go(:,surviving_neurons,:);
n_neurons = sum(surviving_neurons);

%% normalization
total_n_bins = 0;

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    total_n_bins = total_n_bins + epochs.n_bins.(epoch);
end

% preallocation
psth_concat = nan(total_n_bins,n_neurons);

% iterate through task epochs
last_idx = 0;
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    
    % iterate through contrasts
    for jj = 1 : epoch_n_contrasts
        idcs = (1 : epochs.n_bins.(epoch)) + last_idx;
        psth_concat(idcs,:) = psths.(epoch)(:,:,jj);
        last_idx = idcs(end);
    end
end

% z-scoring
mus = nanmean(psth_concat,1);
sigs = nanstd(psth_concat,0,1);
zpsth_concat = (psth_concat - mus) ./ sigs;

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    
    % z-scoring
    zpsths.(epoch) = (psths.(epoch) - mus) ./ sigs;
end

%% plot overall modulation

% figure initialization
fig = figure(figopt,...
    'position',[285,150,1.5392e+03,400],...
    'name',['overall_modulation_',contrast_str]);

% axes initialization
sps = gobjects(n_epochs,1);
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    sps(ii) = subplot(1,n_epochs,ii);
    xlabel(sps(ii),sprintf('%s (ms)',epochs.label.(epoch)));
    set(sps(ii),...
        axesopt.default,...
        'xlim',epochs.roi.(epoch),...
        'xtick',unique([0,epochs.roi.(epoch)]),...
        'ylimspec','tight',...
        'plotboxaspectratiomode','auto');
end
xxtick = unique([0;t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
set(sps([2,4]),...
    'xtick',xxtick,...
    'xticklabel',xxticklabel);
set(sps(2:end),...
    'ycolor','none');
ylabel(sps(1),'Projection onto LD');

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    
    % contrast adjustments
    epoch_contrast_str = epochs.contrast.(epoch);
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    epoch_contrast_mode_idx = find(epoch_contrast_set == mode(epoch_contrasts));
    epoch_contrast_clrs = eval([epoch_contrast_str,'_clrs']);

    % graphical object preallocation
    p = gobjects(epoch_n_contrasts,1);
    
    % iterate through contrasts
    for jj = epoch_n_contrasts : -1 : 1
    
        % project onto linear discriminant
        score = psths.(epoch)(:,:,jj) * lda_mdl.Coeffs(2,1).Linear;

        % plot projection onto linear discriminant
        p(jj) = plot(sps(ii),epochs.time.(epoch),score,...
            'color',epoch_contrast_clrs(jj,:),...
            'linewidth',2);
    end

    % legend
    leg_str = arrayfun(@(x,y,z)sprintf('%s_%s = %i',x,y,z),...
        repmat(upper(epoch_contrast_str(1)),epoch_n_contrasts,1),...
        repmat(epoch_contrast_str(2),epoch_n_contrasts,1),epoch_contrast_set,...
        'uniformoutput',false);
    leg = legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
        'location','north',...
        'box','on');
end

% axes linkage
linkaxes(sps,'y');

% highlight spike integration window used in LDA
lda_roi = post_s2_delay - [spkintegration_window, 0];
xpatch = [lda_roi,fliplr(lda_roi)];
ypatch = [[1,1]*min(ylim(sps(5))),[1,1]*max(ylim(sps(5)))];
p_roi = patch(sps(5),xpatch,ypatch,'k',...
    'facealpha',.15,...
    'edgecolor','none',...
    'handlevisibility','off');
uistack(p_roi,'bottom');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end