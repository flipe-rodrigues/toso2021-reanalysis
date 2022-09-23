%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% plot overall modulation

% figure initialization
fig = figure(figopt,...
    'position',[1,285,1.5392e+03,285],...
    'name',['overall_modulation_',epoch_contrast_str]);

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
epoch_labels = fieldnames(epochs.label);
n_epochs = numel(epoch_labels);
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    
    %
    epoch_contrast_str = epochs.contrast.(epoch);
    epoch_contrasts = eval(epoch_contrast_str);
    epoch_contrast_set = eval([epoch_contrast_str,'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    epoch_contrast_mode_idx = find(epoch_contrast_set == mode(epoch_contrasts));
    epoch_contrast_clrs = eval([epoch_contrast_str,'_clrs']);
%     epoch_contrast_units = eval([epoch_contrast_str(1),'_units']);
%     epoch_contrast_lbl = [upper(epoch_contrast_str(1)),'_',epoch_contrast_str(2)];

    % graphical object preallocation
    p = gobjects(epoch_n_contrasts,1);
    
    % iterate through contrasts
    for jj = epoch_n_contrasts : -1 : 1
    
        % compute modulation stats
        nan_flags = all(isnan(zpsths.(epoch)(:,:,jj)),2);
        mod_mu = nanmean(zpsths.(epoch)(~nan_flags,:,jj),2);
        mod_sig = nanstd(zpsths.(epoch)(~nan_flags,:,jj),0,2);
        mod_sem = mod_sig / sqrt(n_neurons);
        
        % patch s.e.m.
        xpatch = [epochs.time.(epoch)(~nan_flags),...
            fliplr(epochs.time.(epoch)(~nan_flags))];
        ypatch = [mod_mu-mod_sem;flipud(mod_mu+mod_sem)];
        patch(sps(ii),xpatch,ypatch,epoch_contrast_clrs(jj,:),...
            'facealpha',.25,...
            'edgecolor','none');
        
        % plot mean
        p(jj) = plot(sps(ii),epochs.time.(epoch)(~nan_flags),mod_mu,...
            'color',epoch_contrast_clrs(jj,:),...
            'linewidth',1.5);
    end
end

% legend
leg_str = arrayfun(@(x,y,z)sprintf('%s_%s = %i',x,y,z),...
    repmat(upper(epoch_contrast_str(1)),epoch_n_contrasts,1),...
    repmat(epoch_contrast_str(2),epoch_n_contrasts,1),epoch_contrast_set,...
    'uniformoutput',false);
leg = legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
    'position',[0.915,0.425,0.07,0.3],...
    'box','on');

% axes linkage
linkaxes(sps,'y');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end