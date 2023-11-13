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
    ...'t1','i1';...
    ...'t1','i2';...
    ...'t2','i1';...
    ...'t2','i2';...
    };
n_all_pairs = size(all_pair_labels,1);

% iterate through pairs
for ii = 1 : n_all_pairs
    
    % pair specification
    v1_lbl = all_pair_labels{ii,1};
    v2_lbl = all_pair_labels{ii,2};
    v1 = eval(v1_lbl);
    v2 = eval(v2_lbl);
    v_pairs = [v1,v2];
    v_pairset = unique(v_pairs(valid_flags,:),'rows');
    n_v_pairs = size(v_pairset,1);
    
    % preallocation
    n_samples = nan(n_v_pairs,1);
    
    % iterate through current pairs
    for jj = 1 : n_v_pairs
        s_flags = all(v_pairs == v_pairset(jj,:),2);
        trial_flags = ...
            valid_flags & ...
            s_flags;
        if sum(trial_flags) == 0
            continue;
        end
        
        % compute average performance for the current pair
        n_samples(jj) = sum(trial_flags);
    end
    
    % nan filtering
    v_pairset = v_pairset(~isnan(n_samples),:);
    n_samples = n_samples(~isnan(n_samples));
    
    % normalization
    p_sampling = n_samples / max(n_samples);
    
    v1_set = unique(v_pairset(:,1));
    v2_set = unique(v_pairset(:,2));
    
    % figure & axes initialization
    fig = figure(figopt,...
        'name',sprintf('sampling_scheme_%s%s',v1_lbl,v2_lbl));
    axes(axesopt.default,...
        'xlim',tfun([v1_set(1),v1_set(end)]) + [-1,1] * .1 * range(tfun(v1_set)),...
        'ylim',tfun([v2_set(1),v2_set(end)]) + [-1,1] * .1 * range(tfun(v2_set)),...
        'xtick',tfun(v1_set),...
        'ytick',tfun(v2_set),...
        'xticklabel',num2cell(v1_set),...
        'yticklabel',num2cell(v2_set),...
        ...'xticklabelrotation',45,...
        ...'yticklabelrotation',45,...
        'xscale','linear',...
        'yscale','linear');
    if contains(v1_lbl,'t')
        v1_units = s_units;
    elseif contains(v1_lbl,'i')
        v1_units = d_units;
    end
    if contains(v2_lbl,'t')
        v2_units = s_units;
    elseif contains(v2_lbl,'i')
        v2_units = d_units;
    end
    xlabel(sprintf('%s (%s)',upper(insertAfter(v1_lbl,v1_lbl(1),'_')),v1_units));
    ylabel(sprintf('%s (%s)',upper(insertAfter(v2_lbl,v2_lbl(1),'_')),v2_units));
    
    % colorbar
    cbar = colorbar(...
        'color','k',...
        'limits',[min(round(p_sampling,2)),max(round(p_sampling,2))],...
        'box',axesopt.default.box,...
        'linewidth',axesopt.default.linewidth,...
        'tickdirection','out',...
        'ticklength',unique(axesopt.default.ticklength),...
        'fontsize',axesopt.default.fontsize,...
        'ticks',unique(round(p_sampling,2)));
    cbar.Label.String = 'Relative sampling frequency';
    cbar.Label.String = 'Relative frequency';
    cbar.Label.Rotation = -90;
    cbar.Label.VerticalAlignment = 'bottom';
    cbar.TickLabels(2:end-1) = {''};
    
    % plot performance
    scatter(...
        tfun(v_pairset(:,1)),...
        tfun(v_pairset(:,2)),...
        250,p_sampling,'s','filled',...
        'markeredgecolor','k',...
        'linewidth',1.5);
    
    % iterate through pairs
    for jj = 1 : n_v_pairs
        if v_pairset(jj,2) > v_pairset(jj,1)
            if ii == 1
                vertical_gain = .865;
            else
                vertical_gain = .92;
            end
        else
            if ii == 1
                vertical_gain = 1.175;
            else
                vertical_gain = 1.095;
            end
        end
        
        % annotate performance
        text(tfun(v_pairset(jj,1)),tfun(v_pairset(jj,2) * vertical_gain),...
            sprintf('%.2f',p_sampling(jj)),...
            'color','k',...
            'fontsize',axesopt.default.fontsize,...
            'horizontalalignment','center',...
            'verticalalignment','middle');
    end
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end

%% plot marginal distributions

task_vars = unique(all_pair_labels);
n_task_vars = numel(task_vars);

for ii = 1 : n_task_vars
    
    % parse current task variable
    v_lbl = task_vars{ii};
    v = eval(v_lbl);
    v_set = unique(v(valid_flags));
    n_v = numel(v_set);
    
    % figure & axes initialization
    fig = figure(figopt,...
        ...'position',[744,630,460,420],...
        'name',sprintf('marginal_%s',v_lbl));
    
    % axes initialization
    yytick = 0 : 10e3 : 100e3;
    yyticklabel = arrayfun(@(x)sprintf('%ik',x),yytick/1e3,...
        'uniformoutput',false);
    axes(axesopt.default,...
        'plotboxaspectratio',[2,1,1],...
        'position',[0.176 0.15 0.5785 0.65],...
        'xlim',tfun([v_set(1),v_set(end)]) + [-1,1] * .1 * range(tfun(v_set)),...
        'xtick',tfun(v_set),...
        'xticklabel',num2cell(v_set),...
        'ylimspec','tight',...
        'ytick',yytick,...
        'yticklabel',yyticklabel,...
        'yaxislocation','left');
    if contains(v_lbl,'1')
        set(gca,...
            'xaxislocation','top',...
            'ydir','reverse');
    end
    if contains(v_lbl,'t')
        v_units = s_units;
    elseif contains(v_lbl,'i')
        v_units = d_units;
    end
    xlabel(sprintf('%s (%s)',upper(insertAfter(v_lbl,v_lbl(1),'_')),v_units));
    ylabel('Trial count');
    
    % compute marginal distribution
    counts = histcounts(v(valid_flags));
    counts = counts(counts ~= 0);
    
    % iterate through task variable set
    for vv = 1 : n_v
        
        % plot marginal distribution
        barwidth = tfun(range(v_set)) * .01;
        xoffset = barwidth * (ss - round(n_subjects / 2));
        edges = tfun(v_set(vv)) + [-1,1]/2 * barwidth;
        xpatch = [edges,fliplr(edges)];
        ypatch = [zeros(1,2),[1,1]*counts(vv)];
        patch(xpatch,ypatch,stim_clrs(2,:),...
            'facealpha',1,...
            'edgecolor','w',...
            'linewidth',1);
    end
        
    % iterate through subjects
    for ss = 1 : n_subjects
        subject_flags = subjects == ss;
        trial_flags = ...
            valid_flags & ...
            subject_flags;
        
        % compute marginal distribution
        counts = histcounts(v(trial_flags));
        counts = counts(counts ~= 0);
        
        % iterate through task variable set
        for vv = 1 : n_v
            
            % plot marginal distribution
            barwidth = tfun(range(v_set)) * .01;
            xoffset = barwidth * (ss - round(n_subjects / 2));
            edges = tfun(v_set(vv)) + [-1,1]/2 * barwidth;
            xpatch = [edges,fliplr(edges)];
            ypatch = [zeros(1,2),[1,1]*counts(vv)];
            patch(xpatch+xoffset,ypatch,stim_clrs(1,:),...
                'facealpha',1,...
                'edgecolor','w',...
                'linewidth',1);
        end
    end
        
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end