%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% bootstrap settings
n_boots = 50;

%% compute Si-aligned stereotypy

% preallocation
fr_correlation = struct(...
    's1',nan(n_neurons_total,1),...
    's2',nan(n_neurons_total,1));
fr_entropy = struct(...
    's1',nan(n_neurons_total,1),...
    's2',nan(n_neurons_total,1));

% entropy function handle
entropyfun = @(p) -nansum(p .* log2(p));

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
    
    % preallocation
    s1_rhos = nan(n_boots,2,2);
    s2_rhos = nan(n_boots,2,2);
    
    % iterate through bootstrap iterations
    for bb = 1 : n_boots
        train_flags = ismember(1:n_trials,...
            randperm(n_trials,floor(n_trials/2)));
        
        % compute correlation coefficients
        s1_rhos(bb,:,:) = corrcoef(...
            nanmean(s1_spkrates(train_flags,:)),...
            nanmean(s1_spkrates(~train_flags,:)));
        s2_rhos(bb,:,:) = corrcoef(...
            nanmean(s2_spkrates(train_flags,:)),...
            nanmean(s2_spkrates(~train_flags,:)));
    end
    
    % compute stereotypy coefficient
    fr_correlation.s1(nn) = nanmean(s1_rhos(:,1,2));
    fr_correlation.s2(nn) = nanmean(s2_rhos(:,1,2));
    
%     % preallocation
%     s1_rhos = nan(n_trials,2,2);
%     s2_rhos = nan(n_trials,2,2);
%     
%     % iterate through trials
%     for tt = 1 : n_trials
% 
%         % compute correlation coefficients
%         nan_flags = isnan(s1_spkrates(tt,:));
%         s1_rhos(tt,:,:) = corrcoef(s1_spkrates(tt,~nan_flags),...
%             nanmean(s1_spkrates((1:n_trials)~=tt,~nan_flags)));
%         nan_flags = isnan(s2_spkrates(tt,:));
%         s2_rhos(tt,:,:) = corrcoef(s2_spkrates(tt,~nan_flags),...
%             nanmean(s2_spkrates((1:n_trials)~=tt,~nan_flags)));
%     end

    % preallocation
    s1_entropy = nan(n_tbins,1);
    s2_entropy = nan(n_tbins,1);
    
    % iterate through time points
    for tt = 1 : n_tbins
        s1_pdf = histcounts(s1_spkrates(:,tt),30,...
            'normalization','pdf');
        s1_pdf(s1_pdf == 0) = nan;
        s1_entropy(tt) = entropyfun(s1_pdf);
        s2_pdf = histcounts(s2_spkrates(:,tt),30,...
            'normalization','pdf');
        s2_pdf(s2_pdf == 0) = nan;
        s2_entropy(tt) = entropyfun(s2_pdf);
    end

    % compute stereotypy coefficient
    fr_entropy.s1(nn) = 1 ./ nanmean(s1_entropy);
    fr_entropy.s2(nn) = 1 ./ nanmean(s2_entropy);
end

%% statistic selection
stat2plot = fr_correlation;

%% plot firing rate of ramping & non-ramping neurons across task epochs

% figure initialization
fig = figure(figopt,...
    'position',[200 200 560 412.5],...
    'name','ramp_fr_stereotypy');

% epoch settings
epochs = fieldnames(stat2plot);
n_epochs = numel(epochs);

% axes initialization
xxoffset = .25;
xxoffsets = [-1,1] * xxoffset;
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
ylabel('Stereotypy coefficient',...
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
yylim = [-1,1];
% yylim = round(quantile(stat2plot.s2,[0,1])/.5)*.5;
yytick = linspace(yylim(1),yylim(2),5);
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
    plot([1,1]*ee,[min(ylim),1],':k',...
        'handlevisibility','off');
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    [~,pval] = kstest2(distro.(epoch){'ramp'},distro.(epoch){'nonramp'});
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