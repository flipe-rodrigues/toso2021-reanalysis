%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% compute Si-aligned firing rate statistics

% preallocation
fr_mu = struct(...
    's1',nan(n_neurons_total,1),...
    's2',nan(n_neurons_total,1));
fr_min = struct(...
    's1',nan(n_neurons,1),...
    's2',nan(n_neurons,1));
fr_range = struct(...
    's1',nan(n_neurons_total,1),...
    's2',nan(n_neurons_total,1));
fr_dynrange = struct(...
    's1',nan(n_neurons_total,1),...
    's2',nan(n_neurons_total,1));
fr_fano = struct(...
    's1',nan(n_neurons_total,1),...
    's2',nan(n_neurons_total,1));

% iterate through neurons
for nn = 1 : n_neurons_total
    progressreport(nn,n_neurons_total,'computing firing statistics');
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
    
    % S1-aligned spike rates
    s1_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags);
    s1_alignment_flags = ...
        padded_time >= s1_alignment & ...
        padded_time < s1_alignment + t1(spike_flags);
    s1_chunk_flags = ...
        padded_time >= s1_alignment & ...
        padded_time < s1_alignment + t_set(end);
    s1_spkrates = spike_rates';
    s1_spkrates(~s1_alignment_flags') = nan;
    s1_spkrates = reshape(...
        s1_spkrates(s1_chunk_flags'),[n_tbins,n_trials])';
    
    % S2-aligned spike rates
    s2_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi;
    s2_alignment_flags = ...
        padded_time >= s2_alignment & ...
        padded_time < s2_alignment + t2(spike_flags);
    s2_chunk_flags = ...
        padded_time >= s2_alignment & ...
        padded_time < s2_alignment + t_set(end);
    s2_spkrates = spike_rates';
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = reshape(...
        s2_spkrates(s2_chunk_flags'),[n_tbins,n_trials])';
    
    % compute average spike rates across trials
    s1_spkrate = mean(s1_spkrates,'omitnan');
    s2_spkrate = mean(s2_spkrates,'omitnan');
    
    % compute observations weights
    s1_weights = sum(~isnan(s1_spkrates));
    s2_weights = sum(~isnan(s1_spkrates));
    s1_weights = s1_weights ./ nansum(s1_weights);
    s2_weights = s2_weights ./ nansum(s2_weights);
    
    % compute firing rate statistics
    fr_mu.s1(nn) = s1_spkrate * s1_weights';
    fr_mu.s2(nn) = s2_spkrate * s2_weights';
    fr_min.s1(nn) = min(s1_spkrate);
    fr_min.s2(nn) = min(s2_spkrate);
    fr_range.s1(nn) = range(s1_spkrate);
    fr_range.s2(nn) = range(s2_spkrate);
    fr_dynrange.s1(nn) = max(s1_spkrate) / min(s1_spkrate);
    fr_dynrange.s2(nn) = max(s2_spkrate) / min(s2_spkrate);
    fr_fano.s1(nn) = var(s1_spkrate) / mean(s1_spkrate);
    fr_fano.s2(nn) = var(s2_spkrate) / mean(s2_spkrate);
end

%% statistic selection
stat2plot = fr_range;

%% plot firing rate of ramping & non-ramping neurons across task epochs

% figure initialization
fig = figure(figopt,...
    'position',[200 200 560 412.5],...
    'name','ramp_fr_range');

% epoch settings
epochs = fieldnames(stat2plot);
n_epochs = numel(epochs);

% horizontal offset between clusters
xxoffset = .25;
xxoffsets = [-1,1] * xxoffset;

% axes initialization
xxtick = unique((1:n_epochs)+[-1;0;1]*xxoffset);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:n_epochs)) = {''};
xxticklabel(ismember(xxtick,1:n_epochs)) = cellfun(@capitalize,epochs,...
    'uniformoutput',false);
axes(axesopt.default,...
    'plotboxaspectratio',[1,2.25,1],...
    'color','none',...
    'xlim',[1,n_epochs]+[-1,1]*xxoffset*3,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'yaxislocation','right',...
    'clipping','off',...
    'layer','bottom');
xlabel('Stimulus period');
ylabel('Firing rate range (Hz)',...
    'rotation',-90,...
    'verticalalignment','bottom');

% preallocation
distro = struct();
counts = struct();
avg = struct();
err = struct();

% choice of average and error functions
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    distro.(epoch) = {...
        stat2plot.(epoch)(cluster_idcs.(epoch){'ramp'});...
        stat2plot.(epoch)(cluster_idcs.(epoch){'nonramp'})};
    avg.(epoch) = [...
        avgfun(stat2plot.(epoch)(cluster_idcs.(epoch){'ramp'}));...
        avgfun(stat2plot.(epoch)(cluster_idcs.(epoch){'nonramp'}))];
    err.(epoch) = [...
        errfun(stat2plot.(epoch)(cluster_idcs.(epoch){'ramp'}));...
        errfun(stat2plot.(epoch)(cluster_idcs.(epoch){'nonramp'}))];
end

% table conversions
distro = struct2table(distro,...
    'rownames',cluster_labels);
avg = struct2table(avg,...
    'rownames',cluster_labels);
err = struct2table(err,...
    'rownames',cluster_labels);

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    plot(ee+xxoffsets,avg.(epoch),...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    for kk = 1 : n_clusters
        errorbar(ee+xxoffsets(kk),avg.(epoch)(kk),...
            err.(epoch)(kk,1),err.(epoch)(kk,2),...
            'color','k',...
            'marker','o',...
            'markersize',7.5,...
            'markeredgecolor','k',...
            'markerfacecolor',ramp_clrs(kk,:),...
            'linewidth',1.5,...
            'capsize',0);
    end
end

% legend
legend({'ramping','non-ramping'},...
    'autoupdate','off',...
    'box','off',...
    'location','southwest');

% update axes
yymax = ceil(max(ylim)/5) * 5;
yylim = [0,yymax];
yytick = linspace(yylim(1),yylim(2),4);
yyticklabel = num2cell(yytick);
yyticklabel(~ismember(yytick,[0,yylim])) = {''};
set(gca,...
    'ylim',yylim + [-1,1] * .05 * range(yylim),...
    'ytick',yytick,...
    'yticklabel',yyticklabel);

% bin settings
edges = linspace(yylim(1),yylim(2),40);

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    counts.(epoch) = {...
        histcounts(stat2plot.(epoch)(cluster_idcs.(epoch){'ramp'}),edges);...
        histcounts(stat2plot.(epoch)(cluster_idcs.(epoch){'nonramp'}),edges)};
end

% table conversion
counts = struct2table(counts,...
    'rownames',cluster_labels);

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    
    % iterate through clusters
    for kk = n_clusters : -1 : 1
        cluster = cluster_labels{kk};
        xx = counts.(epoch){cluster} / nansum(counts.(epoch){cluster});
        xx = xx / max(xx) * xxoffset * 2 * (-1)^(~iseven(kk)) + ee;
        xx = xx .* [1;1];
        xx = [ee; xx(:); ee];
        yy = edges .* [1;1];
        xpatch = [[1;1]*ee;xx(:)];
        ypatch = [edges([end,1])';yy(:)];
        p = patch(xpatch,ypatch,ramp_clrs(kk,:),...
            'edgecolor','none',...
            'facealpha',1,...
            'linewidth',1.5);
        uistack(p,'bottom');
    end
end

% zero line
plot(xlim,[0,0],':k');

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    xx = [-1,1] * .5 / 3 + ee;
    yy = [1,1] * max(yylim);
    plot([1,1]*ee,[min(ylim),yymax],':k',...
        'handlevisibility','off');
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    [~,pval] = kstest2(distro.(epoch){'ramp'},distro.(epoch){'nonramp'});
    pval = pval * n_epochs;
    if pval < .01
        test_str = '**';
        font_size = 16;
    elseif pval < .05
        test_str = '*';
        font_size = 16;
    else
        test_str = 'n.s.';
        font_size = axesopt.default.fontsize;
    end
    text(mean(xx),mean(yy)-.025*range(ylim),test_str,...
        'color','k',...
        'fontsize',font_size,...
        'horizontalalignment','center',...
        'verticalalignment','bottom');
end

% compare the two stimulus epochs
[~,pval] = kstest2(...
    vertcat(distro.s1{'ramp'},distro.(epoch){'nonramp'}),...
    vertcat(distro.s2{'ramp'},distro.(epoch){'nonramp'}));
fprintf('Two-sample KS p-value (S1 vs. S2 firing rate range): %.2f\n',pval);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end