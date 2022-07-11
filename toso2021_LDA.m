%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% plot overall modulation

% figure initialization
fig = figure(figopt,...
    'position',[1,285,1.5392e+03,500],...
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
ylabel(sps(1),'Firing rate (z-scored)');

% iterate through task epochs
for ii = 1 : n_epochs
    epoch = epoch_labels{ii};
    
    % parse contrast
    epoch_contrast_str = epochs.contrast.(epoch);
    epoch_contrasts = eval(epoch_contrast_str);
    if contains(epoch_contrast_str,'prev')
        epoch_contrast_str = strrep(epoch_contrast_str,'prev_','');
    end
    epoch_contrast_set = eval([epoch_contrast_str(1:end-1),'_set']);
    epoch_n_contrasts = numel(epoch_contrast_set);
    epoch_contrast_clrs = eval([epoch_contrast_str,'_clrs']);
    
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
    
    % legend
    leg_str = arrayfun(@(x,y,z)sprintf('%s_%s = %i',x,y,z),...
        repmat(upper(epoch_contrast_str(1)),epoch_n_contrasts,1),...
        repmat(epoch_contrast_str(2),epoch_n_contrasts,1),epoch_contrast_set,...
        'uniformoutput',false);
    leg = legend(sps(ii),...
        p(isgraphics(p)),leg_str(isgraphics(p)),...
        'box','on');
end

% axes linkage
linkaxes(sps,'y');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%%
