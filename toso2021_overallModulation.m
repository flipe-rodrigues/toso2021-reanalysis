%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% ROI settings

% preallocation
rois = struct();
n_bins = struct();
time = struct();
psths = struct();
zpsths = struct();

% roi definition
rois.pre_s1 = [-500,unique(pre_t1_delay(valid_flags))];
rois.s1 = [0,t_set(end)];
rois.inter_s1s2 = [0,inter_t1t2_delay];
rois.s2 = [0,t_set(end)];
rois.post_s2 = [0,1] * post_t2_delay;
rois.go = [0,450];

% roi alignment labels
alignments.pre_s1 = 'pre-T_{1}-delay';
alignments.s1 = 'T_{1}';
alignments.inter_s1s2 = 'inter T_{1}-T_{2} delay';
alignments.s2 = 'T_{2}';
alignments.post_s2 = 'post-T_{2}-delay';
alignments.go = 'go cue';

% iterate through task epochs
task_epochs = fieldnames(rois);
n_epochs = numel(task_epochs);
for ii = 1 : n_epochs
    epoch = task_epochs{ii};
    n_bins.(epoch) = range(rois.(epoch)) * psthbin;
    time.(epoch) = linspace(rois.(epoch)(1),rois.(epoch)(2),n_bins.(epoch));
    psths.(epoch) = nan(n_bins.(epoch),n_neurons,n_contrasts);
end

%% contrast settings
contrast_str = 't2';
contrasts = eval(contrast_str);
contrast_set = eval([contrast_str(1:end-1),'_set']);
n_contrasts = numel(contrast_set);
contrast_mode_idx = find(contrast_set == mode(contrasts));
contrast_clrs = eval([contrast_str,'_clrs']);

%% construct s2-aligned psths

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
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
            valid_time >= pre_s1_alignment + rois.pre_s1(1) & ...
            valid_time < pre_s1_alignment + pre_t1_delay(spike_flags);
        pre_s1_chunk_flags = ...
            valid_time >= pre_s1_alignment + rois.pre_s1(1) & ...
            valid_time < pre_s1_alignment + rois.pre_s1(2);
        pre_s1_spkrates = spike_rates;
        pre_s1_spkrates(~pre_s1_alignment_flags') = nan;
        pre_s1_spkrates = reshape(...
            pre_s1_spkrates(pre_s1_chunk_flags'),[n_bins.pre_s1,n_trials])';
        
        % T1-aligned spike rates
        s1_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags);
        s1_alignment_flags = ...
            valid_time >= s1_alignment + rois.s1(1) & ...
            valid_time < s1_alignment + t1(spike_flags);
        s1_chunk_flags = ...
            valid_time >= s1_alignment + rois.s1(1) & ...
            valid_time < s1_alignment + rois.s1(2);
        s1_spkrates = spike_rates;
        s1_spkrates(~s1_alignment_flags') = nan;
        s1_spkrates = reshape(...
            s1_spkrates(s1_chunk_flags'),[n_bins.s1,n_trials])';
        
        % inter T1-T2 delay-aligned spike rates
        inter_s1s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags);
        inter_s1s2_alignment_flags = ...
            valid_time >= inter_s1s2_alignment + rois.inter_s1s2(1) & ...
            valid_time < inter_s1s2_alignment + inter_t1t2_delay;
        inter_s1s2_chunk_flags = ...
            valid_time >= inter_s1s2_alignment + rois.inter_s1s2(1) & ...
            valid_time < inter_s1s2_alignment + rois.inter_s1s2(2);
        inter_s1s2_spkrates = spike_rates;
        inter_s1s2_spkrates(~inter_s1s2_alignment_flags') = nan;
        inter_s1s2_spkrates = reshape(...
            inter_s1s2_spkrates(inter_s1s2_chunk_flags'),[n_bins.inter_s1s2,n_trials])';
        
        % T2-aligned spike rates
        s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            inter_t1t2_delay;
        s2_alignment_flags = ...
            valid_time >= s2_alignment + rois.s2(1) & ...
            valid_time < s2_alignment + t2(spike_flags);
        s2_chunk_flags = ...
            valid_time >= s2_alignment + rois.s2(1) & ...
            valid_time < s2_alignment + rois.s2(2);
        s2_spkrates = spike_rates;
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(...
            s2_spkrates(s2_chunk_flags'),[n_bins.s2,n_trials])';
        
        % post T2 delay-aligned spike rates
        post_s2_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            inter_t1t2_delay + ...
            t2(spike_flags);
        post_s2_alignment_flags = ...
            valid_time >= post_s2_alignment + rois.post_s2(1) & ...
            valid_time < post_s2_alignment + post_t2_delay;
        post_s2_chunk_flags = ...
            valid_time >= post_s2_alignment + rois.post_s2(1) & ...
            valid_time < post_s2_alignment + rois.post_s2(2);
        post_s2_spkrates = spike_rates;
        post_s2_spkrates(~post_s2_alignment_flags') = nan;
        post_s2_spkrates = reshape(...
            post_s2_spkrates(post_s2_chunk_flags'),[n_bins.post_s2,n_trials])';
        
        % go-aligned spike rates
        go_alignment = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            inter_t1t2_delay + ...
            t2(spike_flags) + ...
            post_t2_delay;
        go_alignment_flags = ...
            valid_time >= go_alignment + rois.go(1) & ...
            valid_time < go_alignment + rois.go(2);
        go_chunk_flags = go_alignment_flags;
        go_spkrates = spike_rates;
        go_spkrates(~go_alignment_flags') = nan;
        go_spkrates = reshape(...
            go_spkrates(go_chunk_flags'),[n_bins.go,n_trials])';

        % compute mean spike density function
        psths.pre_s1(:,nn,ii) = nanmean(pre_s1_spkrates,1);
        psths.s1(:,nn,ii) = nanmean(s1_spkrates,1);
        psths.inter_s1s2(:,nn,ii) = nanmean(inter_s1s2_spkrates,1);
        psths.s2(:,nn,ii) = nanmean(s2_spkrates,1);
        psths.post_s2(:,nn,ii) = nanmean(post_s2_spkrates,1);
        psths.go(:,nn,ii) = nanmean(go_spkrates,1);
    end
end

%% normalization
total_n_bins = 0;

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = task_epochs{ii};
    total_n_bins = total_n_bins + n_bins.(epoch);
end

% preallocation
psth_concat = nan(total_n_bins,n_neurons);

% iterate through task epochs
last_idx = 0;
for ii = 1 : n_epochs
    epoch = task_epochs{ii};
    
    % iterate through contrasts
    for jj = 1 : n_contrasts
        idcs = (1 : n_bins.(epoch)) + last_idx;
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
    epoch = task_epochs{ii};
    
    % z-scoring
    zpsths.(epoch) = (psths.(epoch) - mus) ./ sigs;
end

%% plot overall modulation

% figure initialization
fig = figure(figopt,...
    'position',[1,285,1.5392e+03,285],...
    'name',['overall_modulation_',contrast_str]);

% axes initialization
sps = gobjects(n_epochs,1);
for ii = 1 : n_epochs
    epoch = task_epochs{ii};
    sps(ii) = subplot(1,n_epochs,ii);
    xlabel(sps(ii),sprintf('%s (ms)',alignments.(epoch)));
    set(sps(ii),...
        axesopt.default,...
        'xlim',rois.(epoch),...
        'xtick',unique([0,rois.(epoch)]),...
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
ylabel(sps(1),'Firing rate (z-scored)');

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = task_epochs{ii};
    
    % iterate through contrasts
    for jj = 1 : n_contrasts
    
        % compute modulation stats
        nan_flags = all(isnan(zpsths.(epoch)(:,:,jj)),2);
        mod_mu = nanmean(zpsths.(epoch)(~nan_flags,:,jj),2);
        mod_sig = nanstd(zpsths.(epoch)(~nan_flags,:,jj),0,2);
        mod_sem = mod_sig / sqrt(n_neurons);
        
        % patch s.e.m.
        xpatch = [time.(epoch)(~nan_flags),fliplr(time.(epoch)(~nan_flags))];
        ypatch = [mod_mu-mod_sem;flipud(mod_mu+mod_sem)];
        patch(sps(ii),xpatch,ypatch,contrast_clrs(jj,:),...
            'facealpha',.25,...
            'edgecolor','none');
        
        % plot mean
        plot(sps(ii),time.(epoch)(~nan_flags),mod_mu,...
            'color',contrast_clrs(jj,:),...
            'linewidth',1.5);
    end
end

% axes linkage
linkaxes(sps,'y');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end