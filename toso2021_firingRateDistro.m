%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct response
roi = [0,t_set(end)];
roi_n_bins = range(roi) * psthbin;

% preallocation
spkrates = nan(n_neurons,n_i);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'fetching spike counts');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % iterate through intensities
    for ii = 1 : n_i
        i2_flags = i2 == i_set(ii);
        
        % flag trials for the current condition
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i2_flags;
        flagged_trials = find(spike_flags);
        if sum(spike_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:)';
        spike_rates = ...
            conv2(kernel.pdf,1,spike_counts,'valid') / psthbin * 1e3;
        n_trials = sum(spike_flags);
        
        % T2-onset-aligned spike rates
        alignment_onset = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi;
        alignment_flags = ...
            valid_time >= alignment_onset + roi(1) & ...
            valid_time < alignment_onset + t2(spike_flags);
        chunk_flags = ...
            valid_time >= alignment_onset + roi(1) & ...
            valid_time < alignment_onset + roi(2);
        nn_spkrates = spike_rates;
        nn_spkrates(~alignment_flags') = nan;
        nn_spkrates = reshape(nn_spkrates(chunk_flags'),[roi_n_bins,n_trials])';
        
        % store average spike rates
        spkrates(nn,ii) = nanmean(nn_spkrates,[1,2]);
    end
end

%% statistics

% anova
[pvalue,~,stats.anova] = anova1(spkrates,...
    num2str(contrast_set),'off');
stats.anova.multcmp = multcompare(stats.anova,'display','off');

%% plot spike rate distribution

% bin settings
n_bins = 21;
binspan = [0,50];
binedges = linspace(binspan(1),binspan(2),n_bins);

% figure initialization
fig = figure(figopt,...
    'name','GLM_firingRateDistro');
axes(axesopt.default,...
    'xlim',binspan+[-1,1]*.05*range(binspan),...
    'ylimspec','tight',...
    'ytick',0,...
    'ycolor','k',...
    'clipping','off',...
    'plotboxaspectratio',[2.25,1,1]);
xlbl = xlabel('Firing rate (Hz)');
ylabel('PDF');

% graphical object preallocation
p = gobjects(n_i,1);

% iterate through intensities
for ii = 1 : n_i
    
    % compute firing rate distribution
    bincounts = histcounts(spkrates(:,ii),...
        'binedges',binedges);
    bincounts = bincounts / nansum(bincounts);
    
    % plot firing rate distribution
    histogram(...
        'binedges',binedges,...
        'bincounts',bincounts,...
        'facealpha',1/n_i,...
        'edgecolor','none',...
        'facecolor',i2_clrs(ii,:),...
        'linewidth',1.5);
    p(ii) = stairs([0,binedges],[0,bincounts,0],...
        'color',i2_clrs(ii,:),...
        'linewidth',1.5);
end

% update y-axis
ylim(ylim+[0,1]*.05*range(ylim));

% iterate through intensities
for ii = 1 : n_i
    
    % compute average firing rate
    mu = nanmean(spkrates(:,ii));
    
    % plot average firing rate
    scatter(mu,max(ylim),85,...
        'linewidth',1.5,...
        'marker','v',...
        'markeredgecolor',i2_clrs(ii,:),....
        'markerfacecolor',i2_clrs(ii,:),...
        'markerfacealpha',1/n_i);
end

% anova outcome
if pvalue < .01
    str = '**';
elseif pvalue < .05
    str = '*';
else
    str = 'n.s.';
end
text(nanmean(spkrates(:)),max(ylim)*1.05,str,...
    'color','k',...
    'fontsize',xlbl.FontSize,...
    'horizontalalignment','center');
% legend
% leg_str = cellfun(@(x,y)sprintf('%s = %i %s',x,y,contrast_units),...
%     repmat({contrast_lbl},n_contrasts,1),num2cell(contrast_set),...
%     'uniformoutput',false);
% legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
%     'box','on');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end