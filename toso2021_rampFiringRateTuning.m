%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% compute temporal tuning across task epochs

% preallocation
fr_tuning = struct();
for ee = 1 : n_cluster_epochs
    epoch = cluster_epochs{ee};
    fr_tuning.(epoch) = nan(n_neurons,1);
end

% iterate through neurons
for nn = 1 : n_neurons_total
    progressreport(nn,n_neurons_total,'parsing neural data');
    neuron_flags = data.NeuronNumb == neuron_idcs(nn);
    spike_flags = ...
        valid_flags & ...
        neuron_flags;
    if sum(spike_flags) == 0 || ~ismember(nn,flagged_neurons)
        continue;
    end
    
    % fetch spike counts & compute spike rates
    spike_rates = data.SDF(spike_flags,:);
    n_trials = size(spike_rates,1);
    
    % S1-onset-aligned spike rates
    s1_onset_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags);
    s1_onset_alignment_flags = ...
        padded_time >= s1_onset_alignment + cluster_roi(1) & ...
        padded_time < s1_onset_alignment + t1(spike_flags);
    s1_onset_chunk_flags = ...
        padded_time >= s1_onset_alignment + cluster_roi(1)& ...
        padded_time < s1_onset_alignment + cluster_roi(2);
    s1_onset_spkrates = spike_rates';
    s1_onset_spkrates(~s1_onset_alignment_flags') = nan;
    s1_onset_spkrates = reshape(...
        s1_onset_spkrates(s1_onset_chunk_flags'),...
        [cluster_roi_n_bins,n_trials])';
    
    % S1-offset-aligned spike rates
    s1_offset_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags);
    s1_offset_alignment_flags = ...
        padded_time >= s1_offset_alignment - t1(spike_flags) & ...
        padded_time < s1_offset_alignment + cluster_roi(2);
    s1_offset_chunk_flags = ...
        padded_time >= s1_offset_alignment + cluster_roi(1)& ...
        padded_time < s1_offset_alignment + cluster_roi(2);
    s1_offset_spkrates = spike_rates';
    s1_offset_spkrates(~s1_offset_alignment_flags') = nan;
    s1_offset_spkrates = reshape(...
        s1_offset_spkrates(s1_offset_chunk_flags'),...
        [cluster_roi_n_bins,n_trials])';
    
    % S2-onset-aligned spike rates
    s2_onset_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi;
    s2_onset_alignment_flags = ...
        padded_time >= s2_onset_alignment + cluster_roi(1) & ...
        padded_time < s2_onset_alignment + t2(spike_flags);
    s2_onset_chunk_flags = ...
        padded_time >= s2_onset_alignment + cluster_roi(1) & ...
        padded_time < s2_onset_alignment + cluster_roi(2);
    s2_onset_spkrates = spike_rates';
    s2_onset_spkrates(~s2_onset_alignment_flags') = nan;
    s2_onset_spkrates = reshape(...
        s2_onset_spkrates(s2_onset_chunk_flags'),...
        [cluster_roi_n_bins,n_trials])';
    
    % S2-offset-aligned spike rates
    s2_offset_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi + ...
        t2(spike_flags);
    s2_offset_alignment_flags = ...
        padded_time >= s2_offset_alignment - t2(spike_flags) & ...
        padded_time < s2_offset_alignment + cluster_roi(2);
    s2_offset_chunk_flags = ...
        padded_time >= s2_offset_alignment + cluster_roi(1) & ...
        padded_time < s2_offset_alignment + cluster_roi(2);
    s2_offset_spkrates = spike_rates';
    s2_offset_spkrates(~s2_offset_alignment_flags') = nan;
    s2_offset_spkrates = reshape(...
        s2_offset_spkrates(s2_offset_chunk_flags'),...
        [cluster_roi_n_bins,n_trials])';
    
    % go-cue-aligned spike rates
    go_cue_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi + ...
        t2(spike_flags) + ...
        post_s2_delay;
    go_cue_alignment_flags = ...
        padded_time >= go_cue_alignment + cluster_roi(1) & ...
        padded_time < go_cue_alignment + cluster_roi(2);
    go_cue_chunk_flags = go_cue_alignment_flags;
    go_cue_spkrates = spike_rates';
    go_cue_spkrates(~go_cue_alignment_flags') = nan;
    go_cue_spkrates = reshape(...
        go_cue_spkrates(go_cue_chunk_flags'),...
        [cluster_roi_n_bins,n_trials])';
       
    % S1-aligned spike rates
    s1_alignment_flags = ...
        padded_time >= s1_onset_alignment & ...
        padded_time < s1_onset_alignment + t1(spike_flags);
    s1_chunk_flags = ...
        padded_time >= s1_onset_alignment & ...
        padded_time < s1_onset_alignment + t_set(end);
    s1_spkrates = spike_rates';
    s1_spkrates(~s1_alignment_flags') = nan;
    s1_spkrates = reshape(...
        s1_spkrates(s1_chunk_flags'),...
        [cluster_roi_n_bins,n_trials])';
     
    % S2-aligned spike rates
    s2_alignment_flags = ...
        padded_time >= s2_onset_alignment & ...
        padded_time < s2_onset_alignment + t2(spike_flags);
    s2_chunk_flags = ...
        padded_time >= s2_onset_alignment & ...
        padded_time < s2_onset_alignment + t_set(end);
    s2_spkrates = spike_rates';
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = reshape(...
        s2_spkrates(s2_chunk_flags'),...
        [cluster_roi_n_bins,n_trials])';
    
    % iterate through cluster epochs
    for ee = 1 : n_cluster_epochs
        epoch = cluster_epochs{ee};
        epoch_spkrates = eval([epoch,'_spkrates']);
        
        % compute average spike rates across trials
        epoch_spkrate = nanmean(epoch_spkrates);
        
        % subtract minimum
        epoch_spkrate = epoch_spkrate - min(epoch_spkrate);
        
        % normalize
        epoch_spkrate = epoch_spkrate ./ nansum(epoch_spkrate);

        % compute center of mass
        fr_tuning.(epoch)(nn) = epoch_spkrate * cluster_roi_time';
    end
end

%% plot temporal tuning of ramping & non-ramping neurons across task epochs

% figure initialization
fig = figure(figopt,...
    'position',[550,350,410,420],...
    'name','ramp_fr_tuning');

% epoch settings
epochs = fieldnames(fr_tuning);
n_epochs = numel(epochs);

% horizontal offset between clusters
xxoffset = .325;
xxoffsets = [-1,1] * xxoffset;

% axes initialization
xxtick = unique((1:n_cluster_epochs)+[-1;0;1]*xxoffset);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:n_cluster_epochs)) = {''};
xxticklabel(ismember(xxtick,1:n_cluster_epochs)) = cellfun(...
    @(x)capitalize(strrep(x,'_',' ')),cluster_epochs,...
    'uniformoutput',false);
yylim = cluster_roi;% / 2;
yytick = linspace(yylim(1),yylim(2),5);
yyticklabel = num2cell(yytick);
yyticklabel(~ismember(yytick,[yylim,0])) = {''};
axes(axesopt.default,...
    'plotboxaspectratio',[2.25,1,1],...
    'color','none',...
    'xlim',[1,n_cluster_epochs]+[-1,1]*xxoffset*2,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'xticklabelrotation',45,...
    'ylim',yylim+[-1,1]*.05*2.25*range(yylim),...
    'ytick',yytick,...
    'yticklabel',yyticklabel,...
    'clipping','off');
xlabel('Task event');
ylabel('Temporal tuning (ms)');

% zero line
plot(xlim,[0,0],':k');

% preallocation
distro = struct();
counts = struct();
counts_up = struct();

% bin settings
edges = linspace(yylim(1),yylim(2),30);

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    distro.(epoch) = {...
        fr_tuning.(epoch)(cluster_idcs.(epoch){'ramp'});...
        fr_tuning.(epoch)(cluster_idcs.(epoch){'nonramp'})};
    counts.(epoch) = {...
        histcounts(fr_tuning.(epoch)(cluster_idcs.(epoch){'ramp'}),edges);...
        histcounts(fr_tuning.(epoch)(cluster_idcs.(epoch){'nonramp'}),edges)};
end

% iterate through alignments
for ee = 1 : n_cluster_epochs
    epoch = cluster_epochs{ee};
    counts_up.(epoch) = {...
        histcounts(fr_tuning.(epoch)(ramp_idcs.(epoch){'up'}),edges)};
end

% table conversions
distro = struct2table(distro,...
    'rownames',cluster_labels);
counts = struct2table(counts,...
    'rownames',cluster_labels);
counts_up = struct2table(counts_up,...
    'rownames',ramp_idcs.Properties.RowNames(1));

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    
    % iterate through clusters
    for kk = n_clusters : -1 : 1
        cluster = cluster_labels{kk};
        xx = counts.(epoch){cluster} / nansum(counts.(epoch){cluster});
        xx = xx / max(xx) * xxoffset * 1.25 * (-1)^(~iseven(kk)) + ee;
        xx = xx .* [1;1];
        xx = [ee; xx(:); ee];
        yy = edges .* [1;1];
        xpatch = [[1;1]*ee;xx(:)];
        ypatch = [edges([end,1])';yy(:)];
        patch(xpatch,ypatch,ramp_clrs(kk,:),...
            'edgecolor','none',...
            'facealpha',1,...
            'linewidth',1.5);
    end
end

% iterate through alignments
for ee = 1 : n_cluster_epochs - 2
    epoch = epochs{ee};
    
    % iterate through clusters
    for kk = 1 % n_clusters : -1 : 1
        xx = counts_up.(epoch){'up'} / nansum(counts.(epoch){'ramp'});
        xx_ref = counts.(epoch){'ramp'} / nansum(counts.(epoch){'ramp'});
        xx = xx / max(xx_ref) * xxoffset * 1.25 * (-1)^(~iseven(kk)) + ee;
        xx = xx .* [1;1];
        xx = [ee; xx(:); ee];
        yy = edges .* [1;1];
        xpatch = [[1;1]*ee;xx(:)];
        ypatch = [edges([end,1])';yy(:)];
        patch(xpatch,ypatch,ramp_clrs(kk,:),...
            'edgecolor','none',...
            'facealpha',1,...
            'linewidth',1.5);
    end
end

% % iterate through alignments
% for ee = 1 : n_epochs
%     epoch = epochs{ee};
%     xx = [-1,1] * .5 / 3 + ee;
%     yy = [1,1] * yylim(2);
%     
%     % iterate through clusters
%     for kk = 1 : n_clusters
%         cluster = cluster_labels{kk};
%         
% %         n_modes = [];
% %         aic = inf;
% %         for mm = 1 : 2
% %             gmfit = fitgmdist(distro.(epoch){cluster},mm);
% %             if aic > gmfit.BIC
% %                 aic = gmfit.BIC;
% %                 n_modes = mm;
% %             end
% %         end
%         
%         [~,pval] = kstest(distro.(epoch){cluster}-mean(distro.(epoch){cluster}));
% %         pval = pval * n_epochs
%         if pval < .01
%             test_str = '**';
%         elseif pval < .05
%             test_str = '*';
%         else
%             test_str = 'n.s.';
%         end
%         text(xx(kk),mean(yy)-.025*range(ylim),test_str,...
%             'color',ramp_clrs(kk,:),...
%             'fontsize',16,...
%             'horizontalalignment','center',...
%             'verticalalignment','bottom');
%     end
% end
% return

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    xx = [-1,1] * .5 / 3 + ee;
    yy = [1,1] * yylim(2);
    plot([1,1]*ee,[min(ylim),yylim(2)],':k',...
        'handlevisibility','off');
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    [~,pval] = kstest2(distro.(epoch){'ramp'},distro.(epoch){'nonramp'});
%     pval = kruskalwallis(vertcat(distro.(epoch){:}),...
%         [ones(size(distro.(epoch){'ramp'}));...
%         zeros(size(distro.(epoch){'nonramp'}))],'off');
    pval = pval * n_epochs;
    if pval < .01
        test_str = '**';
    elseif pval < .05
        test_str = '*';
    else
        test_str = 'n.s.';
    end
    text(mean(xx),mean(yy)-.025*range(ylim),test_str,...
        'color','k',...
        'fontsize',16,...
        'horizontalalignment','center',...
        'verticalalignment','bottom');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end