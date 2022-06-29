%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% sampling scheme

% transfer function
tfun = @(x) log(x);

% pair specification
all_pair_labels = {...
    't1','t2';...
    'i1','i2';...
    't1','i1';...
    't1','i2';...
    't2','i1';...
    't2','i2';...
    };
n_all_pairs = size(all_pair_labels,1);

% iterate through pairs
for ii = 1 : n_all_pairs
    
    % pair specification
    v1_lbl = all_pair_labels{ii,1};
    v2_lbl = all_pair_labels{ii,2};
    v1 = eval(v1_lbl);
    v2 = eval(v2_lbl);
    s_pairs = [v1,v2];
    s_pairset = unique(s_pairs(valid_flags,:),'rows');
    n_pairs = size(s_pairset,1);
    
    % preallocation
    p = nan(n_pairs,1);
    
    % iterate through current pairs
    for jj = 1 : n_pairs
        s_flags = all(s_pairs == s_pairset(jj,:),2);
        trial_flags = ...
            valid_flags & ...
            s_flags;
        if sum(trial_flags) == 0
            continue;
        end
        
        % compute average performance for the current pair
        p(jj) = sum(trial_flags);
    end
    
    % nan filtering
    s_pairset = s_pairset(~isnan(p),:);
    p = p(~isnan(p));
    
    % normalization
    p = p / min(p);
    
    v1_set = unique(s_pairset(:,1));
    v2_set = unique(s_pairset(:,2));
    
    % figure & axes initialization
    fig = figure(figopt,...
        'name','sampling_scheme');
    axes(axesopt.default,...
        'xlim',tfun([v1_set(1),v1_set(end)]) + [-1,1] * .1 * range(tfun(v1_set)),...
        'ylim',tfun([v2_set(1),v2_set(end)]) + [-1,1] * .1 * range(tfun(v2_set)),...
        'xtick',tfun(v1_set),...
        'ytick',tfun(v2_set),...
        'xticklabel',num2cell(v1_set),...
        'yticklabel',num2cell(v2_set),...
        'xticklabelrotation',45,...
        'yticklabelrotation',45,...
        'xscale','linear',...
        'yscale','linear');
    if contains(v1_lbl,'t')
        v1_units = t_units;
    elseif contains(v1_lbl,'i')
        v1_units = i_units;
    end
    if contains(v2_lbl,'t')
        v2_units = t_units;
    elseif contains(v2_lbl,'i')
        v2_units = i_units;
    end
    xlabel(sprintf('%s (%s)',upper(v1_lbl),v1_units));
    ylabel(sprintf('%s (%s)',upper(v2_lbl),v2_units));
    
    % colorbar
    cbar = colorbar(...
        'limits',[min(round(p,2)),max(round(p,2))],...
        'box',axesopt.default.box,...
        'linewidth',axesopt.default.linewidth,...
        'tickdirection','out',...
        'ticklength',unique(axesopt.default.ticklength),...
        'fontsize',axesopt.default.fontsize,...
        'ticks',unique(round(p,2)));
    cbar.Label.String = 'Min-normalized frequency';
    cbar.Label.Rotation = -90;
    cbar.Label.VerticalAlignment = 'bottom';
    
    % plot reference lines
    plot(xlim,ylim,':k',...
        'linewidth',1);
    plot(tfun([1,1]*median(s1(valid_flags))),ylim,':k',...
        'linewidth',1);
    plot(xlim,tfun([1,1]*median(s2(valid_flags))),':k',...
        'linewidth',1);
    
    % plot performance
    scatter(...
        tfun(s_pairset(:,1)),...
        tfun(s_pairset(:,2)),...
        250,p,'s','filled',...
        'markeredgecolor','k',...
        'linewidth',1.5);
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end