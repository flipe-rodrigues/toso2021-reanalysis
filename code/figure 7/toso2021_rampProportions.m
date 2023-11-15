%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% plot proportion of ramping & non-ramping neurons across task epochs

% figure initialization
fig = figure(figopt,...
    'position',[550,350,380,420],...
    'name','ramp_proportions');

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
yytick = linspace(0,1,5);
yyticklabel = num2cell(yytick);
yyticklabel(~ismember(yytick,[0,1])) = {''};
axes(axesopt.default,...
    'plotboxaspectratio',[2.25,1,1],...
    'color','none',...
    'xlim',[1,n_cluster_epochs]+[-1,1]*xxoffset*2,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'xticklabelrotation',45,...
    'ylim',[0,1]+[-1,1]*.05*2.25,...
    'ytick',yytick,...
    'yticklabel',yyticklabel,...
    'clipping','off');
xlabel('Task event');
ylabel('Proportion of neurons');

% preallocation
proportion = struct();
proportion_ud = struct();

% iterate through alignments
for ee = 1 : n_cluster_epochs
    epoch = cluster_epochs{ee};
    proportion.(epoch) = [...
        numel(cluster_idcs.(epoch){'ramp'}); numel(cluster_idcs.(epoch){'nonramp'})] ./ ...
        (numel(cluster_idcs.(epoch){'ramp'}) + numel(cluster_idcs.(epoch){'nonramp'}));
end

% iterate through alignments
for ee = 1 : n_cluster_epochs
    epoch = cluster_epochs{ee};
    proportion_ud.(epoch) = [...
        numel(ramp_idcs.(epoch){'up'}); numel(ramp_idcs.(epoch){'down'})] ./ ...
        (numel(cluster_idcs.(epoch){'ramp'}) + numel(cluster_idcs.(epoch){'nonramp'}));
end

% table conversions
proportion = struct2table(proportion,...
    'rownames',cluster_labels);
proportion_ud = struct2table(proportion_ud,...
    'rownames',ramp_idcs.Properties.RowNames(1:2));

% graphics object preallocation
p = gobjects(3,1);

% plot reference lines
plot(xlim,[1,1]*.5,':k',...
    'handlevisibility','off');

% iterate through alignments
for ee = 1 : n_cluster_epochs
    epoch = cluster_epochs{ee};
    plot([1,1]*ee,ylim,':k',...
        'handlevisibility','off');
    xxrange = xxtick(ee*3+(-2:0));
    for kk = 1 : n_clusters
        xx = xxrange(kk*2-1) + [0,1] * xxoffset * (-1)^(kk == 2) * .75;
        xpatch = [xx,fliplr(xx)];
        ypatch = sort(repmat([0,proportion.(epoch)(kk)],1,2));
        p(kk) = patch(xpatch,ypatch,ramp_clrs(kk,:),...
            'facealpha',1,...
            'edgecolor',ramp_clrs(kk,:),...
            'linewidth',1.5,...
            'handlevisibility','off');
    end
end

% iterate through alignments
for ee = 1 : n_cluster_epochs - 2
    epoch = cluster_epochs{ee};
    plot([1,1]*ee,ylim,':k',...
        'handlevisibility','off');
    xxrange = xxtick(ee*3+(-2:0));
    for kk = 1 %: n_clusters
        xx = xxrange(kk*2-1) + [0,1] * xxoffset * (-1)^(kk == 2) * .75;
        xpatch = [xx,fliplr(xx)];
        ypatch = sort(repmat([0,proportion_ud.(epoch)('down')],1,2));
        p(n_clusters + 1) = patch(xpatch,ypatch,'w',...
            'facealpha',1,...
            'edgecolor',ramp_clrs(kk,:),...
            'linewidth',1.5);
    end
end

% legend
legend(p,{'down-ramping / ramping','non-ramping','up-ramping'},...
    'autoupdate','off',...
    'box','off',...
    'location','best');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end