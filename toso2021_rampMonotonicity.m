%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% compute Si-aligned monotonicity

% monotonicity definition
monotonicityfun = @(x) abs(sum(sign(diff(x)) ./ (length(x) - 1)));

% preallocation
monotonicity = struct(...
    's1',nan(n_neurons_total,1),...
    's2',nan(n_neurons_total,1));

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

    % compute firing rate statistics
    monotonicity.s1(nn) = monotonicityfun(s1_spkrate);
    monotonicity.s2(nn) = monotonicityfun(s2_spkrate);
end

%% plot monotonicity of ramping & non-ramping neurons across task epochs

% figure initialization
fig = figure(figopt,...
    'position',[200 200 560 412.5],...
    'name','ramp_monotonicity');

% epoch settings
epochs = fieldnames(monotonicity);
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
    'clipping','off',...
    'layer','bottom');
xlabel('Stimulus period');
ylabel('Monotonicity');

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
        monotonicity.(epoch)(cluster_idcs.(epoch){'ramp'});...
        monotonicity.(epoch)(cluster_idcs.(epoch){'nonramp'})};
    avg.(epoch) = [...
        avgfun(monotonicity.(epoch)(cluster_idcs.(epoch){'ramp'}));...
        avgfun(monotonicity.(epoch)(cluster_idcs.(epoch){'nonramp'}))];
    err.(epoch) = [...
        errfun(monotonicity.(epoch)(cluster_idcs.(epoch){'ramp'}));...
        errfun(monotonicity.(epoch)(cluster_idcs.(epoch){'nonramp'}))];
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
yymax = 1;ceil(max(ylim)/5) * 5;
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
        histcounts(monotonicity.(epoch)(cluster_idcs.(epoch){'ramp'}),edges);...
        histcounts(monotonicity.(epoch)(cluster_idcs.(epoch){'nonramp'}),edges)};
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
    %     pval = kruskalwallis(vertcat(distro.(epoch){:}),[],'off');
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