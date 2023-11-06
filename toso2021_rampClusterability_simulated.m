%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% gap statistic settings
max_k = 15;
min_pcs2consider = 2;
max_pcs2consider = 2;
n_pcs2consider = max_pcs2consider - min_pcs2consider + 1;
pc_clrs = flipud(gray(max_pcs2consider+1));

%% temporal smoothing kernel
gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);

%% color scheme
simramp_clrs = [[0,1,1]*.85; ramp_clrs(2,:)];

%% ROI settings
t_units = 1e3;
ti = 0;
tf = t_set(end);
roi = [ti,tf];
T = (tf - ti) / psthbin;
t = linspace(ti,tf,T);
dur = (tf - ti) / t_units;
ti_padded = ti + gauss_kernel.paddx(1);
tf_padded = tf + gauss_kernel.paddx(2);
T_padded = (tf_padded - ti_padded) / psthbin;
t_padded = linspace(ti_padded,tf_padded,T_padded);
dur_padded = (tf_padded - ti_padded) / t_units;
dt = diff(t(1:2));

%% simulation settings
N = n_neurons;
K = 100;
M = 1;
ramp_proportion = proportion.s2(1);

%% firing rate statistics
fr.min = fr_min.s2(flagged_neurons);
fr.range = fr_range.s2(flagged_neurons);

%% simulate spike rates

% preallocation
R = nan(T,N,K);
L = categorical(zeros(N,1),[0,1],{'ramp','nonramp'});
L2 = categorical(zeros(N,1),[0,1,2],{'up','non','down'});

% iterate through neurons
for nn = 1 : N
    progressreport(nn,N,'simulating')
    
    % label assignment
    if rand > ramp_proportion
        L(nn) = 'nonramp';
        L2(nn) = 'non';
        m = randi(M);
        mus = unifrnd(ti,tf,m,1);
        sigmas = unifrnd(.15,.3,m,1) * (tf - ti);
        x_padded = sum(normpdf(t_padded,mus,sigmas),1) * randn;
    else
        L(nn) = 'ramp';
        slope = (-1)^(rand>.5) + randn;
        if slope > 0
            L2(nn) = 'up';
        else
            L2(nn) = 'down';
        end
        x_padded = t_padded .* slope;
    end
    
    % match firing statistics to those found in striatal data
    x_padded = normalize01(x_padded,2) * fr.range(nn) + fr.min(nn);
    
    % iterate through trials
    for kk = 1 : K
        [~,ts] = poissonprocess(x_padded,dur_padded);
        spk_times = ts * t_units + ti_padded;
        spk_counts = histcounts(spk_times,'binedges',t_padded);
        R(:,nn,kk) = conv(spk_counts/(dt/t_units),gauss_kernel.pdf,'valid');
    end
end

% compute cross-trial averages
r = nanmean(R,3);

%% parse cluster indices & colors
ramp_flags = L == 'ramp';
nonramp_flags = L == 'nonramp';
upramp_flags = L2 == 'up';
downramp_flags = L2 == 'down';
cluster_ids = 2 - ramp_flags';
cluster_clrs = simramp_clrs(cluster_ids,:);

%% nanify to mimic trials ending at different timepoints ???


%% normalization
mus = nanmean(r,[1,3]);
sigs = nanstd(r,0,[1,3]);
z = (r - mus) ./ sigs;

%% PCA
pca_neuron_coeff = pca(zscore(r,0,1));
[~,pca_time_score] = pca(zscore(r,0,2)');
[thetas,~] = cart2pol(pca_neuron_coeff(:,1),pca_neuron_coeff(:,2));

%% t-SNE
tsne_embeddings = tsne(zscore(r,0,2)');

%% compute pairwise distances between PC coefficients

% all to all pairwise distances
dissimilarity_mat = pdist2(...
    pca_neuron_coeff(:,1:2),pca_neuron_coeff(:,1:2));
triu_mask = dissimilarity_mat == triu(dissimilarity_mat);
dissimilarity = mat2colvec(dissimilarity_mat(triu_mask));

% ramp to all pairwise distances
dissimilarity_mat_ramp = pdist2(...
    pca_neuron_coeff(ramp_flags,1:2),pca_neuron_coeff(ramp_flags,1:2));
triu_mask_ramp = dissimilarity_mat_ramp == triu(dissimilarity_mat_ramp);
dissimilarity_ramp = mat2colvec(dissimilarity_mat_ramp(triu_mask_ramp));

% non-ramp to all pairwise distances
dissimilarity_mat_non = pdist2(...
    pca_neuron_coeff(nonramp_flags,1:2),pca_neuron_coeff(nonramp_flags,1:2));
triu_mask_non = dissimilarity_mat_non == triu(dissimilarity_mat_non);
dissimilarity_non = mat2colvec(dissimilarity_mat_non(triu_mask_non));

%% compute monotonicity
monotonicityfun = @(x) sum(sign(diff(x)) ./ (size(x,1) - 1));
monotonicity = monotonicityfun(r);

%% linear regression

% preallocation
slopes = nan(N,1);
r2 = nan(N,1);

% iterate through neurons
for nn = 1 : N
    mdl = fitlm(x,z(:,nn));
    slopes(nn) = mdl.Coefficients.Estimate(2);
    r2(nn) = mdl.Rsquared.Ordinary;
end

%% clusterability metrics structure

% preallocation
clusterability = struct();

% assignment
% clusterability.dissimilarity = dissimilarity;
clusterability.thetas = thetas;
clusterability.slopes = slopes;

% parse clusterability metrics
clusterability_metrics = fieldnames(clusterability);
n_metrics = numel(clusterability_metrics);

%% hartigan's dip test

% preallocation
diptestpval = struct();

% iterate through clusterability metrics
for mm = 1 : n_metrics
    metric = clusterability_metrics{mm};
    
    % store dip test p-values
    [~,diptestpval.(metric)] = diptest(clusterability.(metric));
end

%% compute gap statistics

% preallocation
cluster_eval = struct();

% iterate through number of PCs considered
for pc = min_pcs2consider : max_pcs2consider
    progressreport(pc-min_pcs2consider+1,n_pcs2consider,...
        'computing gap statistics');
    
    % evaluate clustering based on neuron-wise PC coefficients
    curr_eval = evalclusters(pca_neuron_coeff(:,1:pc),'kmeans','gap',...
        'referencedistribution','uniform',...
        'klist',1:max_k);
    cluster_eval.neuron.gap(pc,:) = curr_eval.CriterionValues;
    cluster_eval.neuron.se(pc,:) = curr_eval.SE;
    cluster_eval.neuron.k(pc) = curr_eval.OptimalK;
    
    % evaluate clustering based on time-wise PC scores
    curr_eval = evalclusters(pca_time_score(:,1:pc),'kmeans','gap',...
        'referencedistribution','uniform',...
        'klist',1:max_k);
    cluster_eval.time.gap(pc,:) = curr_eval.CriterionValues;
    cluster_eval.time.se(pc,:) = curr_eval.SE;
    cluster_eval.time.k(pc) = curr_eval.OptimalK;
end

%% tiling

% figure initialization
fig = figure(figopt,...
    'name','tiling_simulated');

% axes initialization
xxtick = unique([0;roi';t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
axes(axesopt.default,...
    'xlim',roi,...
    'ylim',[1,N],...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ytick',[1,N],...
    'colormap',hot(2^8),...
    'clipping','off');
title('Simulated PSTH raster');
xlabel('Time (ms) 1_1');
ylabel({'Simulated neuron #','(sorted by PCs)'});

% sort by angular position in PC space
[~,theta_idcs] = sortrows(thetas);
% theta_idcs = (circshift(theta_idcs,75));
sorted_idcs = (circshift(theta_idcs,115));

% color limits
clim = [-2,4];

% plot psth raster
% downramp_idcs = find(downramp_flags);
% nonramp_idcs = find(nonramp_flags);
% upramp_idcs = find(upramp_flags);
% sorted_idcs = [downramp_idcs;nonramp_idcs(theta_idcs);upramp_idcs];
imagesc(roi,[1,N],z(:,sorted_idcs)',clim);

% plot cluster affordances
plot(ti-.1*range(xlim),find(nonramp_flags(sorted_idcs)),...
    'color',simramp_clrs(2,:),...
    'marker','.',...
    'markersize',5);
plot(ti-.085*range(xlim),find(ramp_flags(sorted_idcs)),...
    'color',simramp_clrs(1,:),...
    'marker','.',...
    'markersize',5);

% iterate through clusters
for kk = n_clusters : -1 : 1
    cluster_flags = L == cluster_labels{kk};
    neuron_edges = linspace(.5,n_neurons+.5,50);
    neuron_pdf = histcounts(find(cluster_flags(sorted_idcs)),neuron_edges);
    xx = neuron_pdf / max(neuron_pdf);
    xx = xx .* [1;1] * range(xlim) * .035 * -1;
    xx = [0; xx(:); 0];
    yy = neuron_edges .* [1;1];
    xpatch = [[1;1]; xx(:)];
    ypatch = [neuron_edges([end,1])'; yy(:)];
    patch(xpatch,ypatch,simramp_clrs(kk,:),...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
end

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = 'Firing rate (z-score)';
clrlabel.fontsize = axesopt.default.fontsize * 1.1;
clrlabel.rotation = 270;
clrlabel.position = [3.5,sum(clim)/2,0];
set(clrbar,...
    axesopt.colorbar,...
    'color','k',...
    'fontsize',axesopt.default.fontsize);
set(clrbar.Label,...
    'color','k',...
    clrlabel);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% PC neuron-wise coefficient scatter
fig = figure(figopt,...
    'name','pc_neuronwiseCoefficients2D_simulated');
axes(axesopt.default,...
    'colormap',1-gray(2^8),...
    'clipping','off');
xlabel('Neuron-wise PC 1 coefficient_{1}');
ylabel('Neuron-wise PC 2 coefficient_{2}');

% background joint distribution
nbins = 15;
xbounds = quantile(pca_neuron_coeff(:,1),[0,1]);
ybounds = quantile(pca_neuron_coeff(:,2),[0,1]);
xedges = linspace(xbounds(1),xbounds(2),nbins+1);
yedges = linspace(ybounds(1),ybounds(2),nbins+1);
rampcounts2 = histcounts2(...
    pca_neuron_coeff(ramp_flags,1),pca_neuron_coeff(ramp_flags,2),...
    'xbinedges',xedges,...
    'ybinedges',yedges);
noncounts2 = histcounts2(...
    pca_neuron_coeff(nonramp_flags,1),pca_neuron_coeff(nonramp_flags,2),...
    'xbinedges',xedges,...
    'ybinedges',yedges);
C = cat(3,rampcounts2',noncounts2');
P = tensor2rgb(C,simramp_clrs,245/255);
imagesc(...
    xbounds+[1,-1]*diff(xedges(1:2))/2,...
    ybounds+[1,-1]*diff(yedges(1:2))/2,P);

% coefficient scatter
% grapeplot(pca_neuron_coeff(:,1),pca_neuron_coeff(:,2),...
%     'markersize',3,...
%     'markeredgecolor',[0,0,0],...
%     'markerfacecolor',[1,1,1]);
grapeplot(pca_neuron_coeff(nonramp_flags,1),pca_neuron_coeff(nonramp_flags,2),...
    'markersize',3,...
    'markeredgecolor',cluster_clrs(nonramp_flags,:),...
    'markerfacecolor',cluster_clrs(nonramp_flags,:)*.85);
grapeplot(pca_neuron_coeff(ramp_flags,1),pca_neuron_coeff(ramp_flags,2),...
    'markersize',3,...
    'markeredgecolor',cluster_clrs(ramp_flags,:),...
    'markerfacecolor',cluster_clrs(ramp_flags,:)*.85);

% update axis
xxlim = xlim;
yylim = ylim;
xlim(xxlim + [-1,1] * .05 * range(xxlim));
ylim(yylim + [-1,1] * .05 * range(yylim));
set(gca,...
    'xtick',unique([xxlim,0]),...
    'ytick',unique([yylim,0]),...
    'xticklabel',{'','0',''},...
    'yticklabel',{'','0',''},...
    'xcolor','k',...
    'ycolor','k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% PC time-wise score scatter
fig = figure(figopt,...
    'name','pc_timewiseScores2D_simulated');
axes(axesopt.default,...
    'clipping','on');
xlabel('Time-wise PC 1 score_{1}');
ylabel('Time-wise PC 2 score_{2}');

% background joint distribution
nbins = 15;
xbounds = quantile(pca_time_score(:,1),[0,1]);
ybounds = quantile(pca_time_score(:,2),[0,1]);
xedges = linspace(xbounds(1),xbounds(2),nbins+1);
yedges = linspace(ybounds(1),ybounds(2),nbins+1);
rampcounts2 = histcounts2(...
    pca_time_score(ramp_flags,1),pca_time_score(ramp_flags,2),...
    'xbinedges',xedges,...
    'ybinedges',yedges);
noncounts2 = histcounts2(...
    pca_time_score(nonramp_flags,1),pca_time_score(nonramp_flags,2),...
    'xbinedges',xedges,...
    'ybinedges',yedges);
C = cat(3,rampcounts2',noncounts2');
P = tensor2rgb(C,simramp_clrs);
imagesc(...
    xbounds+[1,-1]*diff(xedges(1:2))/2,...
    ybounds+[1,-1]*diff(yedges(1:2))/2,P);

% coefficient scatter
grapeplot(pca_time_score(nonramp_flags,1),pca_time_score(nonramp_flags,2),...
    'markersize',3,...
    'markeredgecolor',cluster_clrs(nonramp_flags,:),...
    'markerfacecolor',cluster_clrs(nonramp_flags,:)*.85);
grapeplot(pca_time_score(ramp_flags,1),pca_time_score(ramp_flags,2),...
    'markersize',3,...
    'markeredgecolor',cluster_clrs(ramp_flags,:),...
    'markerfacecolor',cluster_clrs(ramp_flags,:)*.85);

% update axis
xxlim = xlim; % quantile(s2_score_time(:,1),[0,1]+[1,-1]*.00);
yylim = ylim; % quantile(s2_score_time(:,2),[0,1]+[1,-1]*.00);
xlim(xxlim + [-1,1] * .05 * range(xxlim));
ylim(yylim + [-1,1] * .05 * range(yylim));
set(gca,...
    'xtick',unique([xxlim,0]),...
    'ytick',unique([yylim,0]),...
    'xticklabel',{'','0',''},...
    'yticklabel',{'','0',''},...
    'xcolor','k',...
    'ycolor','k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% t-SNE embeddings
fig = figure(figopt,...
    'name','tsne_simulated');
axes(axesopt.default,...
    'xscale','linear',...
    'yscale','linear',...
    'clipping','off');
xlabel('t-SNE 1_{1}');
ylabel('t-SNE 2_{2}');

% coefficient scatter
grapeplot(tsne_embeddings(nonramp_flags,1),tsne_embeddings(nonramp_flags,2),...
    'markeredgecolor',cluster_clrs(nonramp_flags,:),...
    'markerfacecolor',cluster_clrs(nonramp_flags,:)*.85);
grapeplot(tsne_embeddings(ramp_flags,1),tsne_embeddings(ramp_flags,2),...
    'markeredgecolor',cluster_clrs(ramp_flags,:),...
    'markerfacecolor',cluster_clrs(ramp_flags,:)*.85);

% update axis
xxlim = xlim;
yylim = ylim;
xlim(xxlim + [-1,1] * .05 * range(xxlim));
ylim(yylim + [-1,1] * .05 * range(yylim));
set(gca,...
    'xtick',unique([xxlim,0]),...
    'ytick',unique([yylim,0]),...
    'xticklabel',{'','0',''},...
    'yticklabel',{'','0',''},...
    'xcolor','k',...
    'ycolor','k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% distributions of clusterability metrics

% figure initialization
fig = figure(figopt,...
    'position',[200,200,560,415],...
    'name','clusterability_simulated');

% axes initialization
n_sps = n_metrics;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(n_sps,1,ii);
    ylabel(sps(ii),'PDF');
end
set(sps,axesopt.default,...
    'plotboxaspectratio',[2.75,1,1],...
    'ticklength',axesopt.default.ticklength,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'ytick',0,...
    'clipping','on');
xlabel(sps(1),'Polar angle of PC coefficients');
xlabel(sps(2),'Linear regression slopes');
% xlabel(sps(3),'Pairwise distance of PC coefficients');

% bin settings
n_bins = 30;

% iterate through clusterability metrics
for mm = 1 : n_metrics
    metric = clusterability_metrics{mm};
    
    % compute distributions of clusterability metric
    bounds.(metric) = quantile(clusterability.(metric),[0,1]);
    edges.(metric) = linspace(bounds.(metric)(1),bounds.(metric)(2),n_bins);
	counts.(metric) = histcounts(clusterability.(metric),edges.(metric));
    
    % plot distribution
    histogram(sps(mm),...
        'binedges',edges.(metric),...
        'bincounts',counts.(metric),...
        'facecolor','w',...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    stairs(sps(mm),edges.(metric),[counts.(metric),0],...
        'color','k',...
        'linewidth',1.5);
    
    % affordances for statistical tests
    pval = diptestpval.(metric);
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
    text(sps(mm),.5,.95,test_str,...
        'color','k',...
        'fontsize',font_size,...
        'horizontalalignment','center',...
        'verticalalignment','bottom',...
        'units','normalized');
end

% update axes
set(sps(1),...
    'xlim',bounds.thetas,...
    'xtick',sort([0,bounds.thetas]),...
    'xticklabel',{'-\pi','0','\pi'});
set(sps(2),...
    'xlim',bounds.slopes,...
    'xtick',sort([0,bounds.slopes]),...
    'xticklabel',num2cell(round(sort([0,bounds.slopes]),2)));
% set(sps(3),...
%     'xlim',bounds.dissimilarity,...
%     'xtick',bounds.dissimilarity,...
%     'xticklabel',num2cell(round(bounds.dissimilarity,2)));

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% distributions of clusterability metrics (split by cluster)

% figure initialization
fig = figure(figopt,...
    'position',[200,200,560,415],...
    'name','clusterability_clustersplit_simulated');

% axes initialization
n_sps = n_metrics;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(n_sps,1,ii);
    ylabel(sps(ii),'PDF');
end
set(sps,axesopt.default,...
    'plotboxaspectratio',[2.75,1,1],...
    'ticklength',axesopt.default.ticklength,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'ytick',0,...
    'clipping','on');
set(sps(1),'ylim',[0,max(counts.thetas)]);
set(sps(2),'ylim',[0,max(counts.slopes)]);
% set(sps(3),'ylim',[0,max(counts.dissimilarity)]);
xlabel(sps(1),'Polar angle of PC coefficients');
xlabel(sps(2),'Linear regression slopes');
% xlabel(sps(3),'Pairwise distance of PC coefficients');

% iterate through clusterability metrics
for mm = 1 : n_metrics
    metric = clusterability_metrics{mm};
    
    % compute distributions of clusterability metric
	counts_ramp = histcounts(clusterability.(metric)(ramp_flags),edges.(metric));
	counts_non = histcounts(clusterability.(metric)(nonramp_flags),edges.(metric));
    
    % plot distribution
    histogram(sps(mm),...
        'binedges',edges.(metric),...
        'bincounts',counts_non,...
        'facecolor',simramp_clrs(2,:),...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    stairs(sps(mm),edges.(metric),[counts_non,0],...
        'color','k',...
        'linewidth',1.5);
    histogram(sps(mm),...
        'binedges',edges.(metric),...
        'bincounts',counts_ramp,...
        'facecolor',simramp_clrs(1,:),...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    histogram(sps(mm),...
        'binedges',edges.(metric),...
        'bincounts',counts_non,...
        'facecolor',simramp_clrs(2,:),...
        'edgecolor','none',...
        'facealpha',.5,...
        'linewidth',1.5);
    stairs(sps(mm),edges.(metric),[counts_ramp,0],...
        'color','k',...
        'linewidth',1.5);
    
    % affordances for statistical tests
    pval = diptestpval.(metric);
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
    text(sps(mm),.5,.95,test_str,...
        'color','k',...
        'fontsize',font_size,...
        'horizontalalignment','center',...
        'verticalalignment','bottom',...
        'units','normalized');
end

% update axes
set(sps(1),...
    'xlim',bounds.thetas,...
    'xtick',sort([0,bounds.thetas]),...
    'xticklabel',{'-\pi','0','\pi'});
set(sps(2),...
    'xlim',bounds.slopes,...
    'xtick',sort([0,bounds.slopes]),...
    'xticklabel',num2cell(round(sort([0,bounds.slopes]),2)));
% set(sps(3),...
%     'xlim',bounds.dissimilarity,...
%     'xtick',bounds.dissimilarity,...
%     'xticklabel',num2cell(round(bounds.dissimilarity,2)));

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot gap statistics
clustering_criteria = fieldnames(cluster_eval);
n_clustering_criteria = numel(clustering_criteria);

% iterate through clustering criteria
for cc = 1 : n_clustering_criteria
    criterion = clustering_criteria{cc};
    
    % figure initialization
    fig = figure(figopt,...
        'name',sprintf('gap_%s_simulated',criterion));
    axes(axesopt.default,...
        'xlim',[1,max_k]+[-1,1]*.05*max_k,...
        'xtick',1:max_k,...
        'ylimspec','tight');
    xlabel('Number of clusters');
    ylabel('Gap value');
    
    % iterate through number of PCs considered
    for pc = min_pcs2consider : max_pcs2consider
        errorbar(1:max_k,...
            cluster_eval.(criterion).gap(pc,:),...
            cluster_eval.(criterion).se(pc,:),...
            'color',pc_clrs(pc,:),...
            'linewidth',1.5,...
            'capsize',0);
        plot(cluster_eval.(criterion).k(pc),...
            cluster_eval.(criterion).gap(pc,cluster_eval.(criterion).k(pc)),...
            'color','k',...
            'marker','o',...
            'markersize',6.5,...
            'markeredgecolor',pc_clrs(pc,:),...
            'markerfacecolor','w',...
            'linewidth',1.5,...
            'handlevisibility','off');
    end
    
    % legend
    leg_str = arrayfun(...
        @(x)sprintf('first %i PCs_{s2}',x),min_pcs2consider:max_pcs2consider,...
        'uniformoutput',false);
    % legend(leg_str);
    
    % update axis
    set(gca,...
        'ytick',ylim,...
        'ylim',ylim+[-1,1]*.05*range(ylim),...
        'yticklabel',{'0',''});
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end

%%
return;
nbins = 15;
xedges = linspace(...
    min(pca_neuron_coeff(:,1),[],'all'),...
    max(pca_neuron_coeff(:,1),[],'all'),nbins+1);
yedges = linspace(...
    min(pca_neuron_coeff(:,2),[],'all'),...
    max(pca_neuron_coeff(:,2),[],'all'),nbins+1);
counts2 = histcounts2(...
    pca_neuron_coeff(:,1),pca_neuron_coeff(:,2),...
    'xbinedges',xedges,...
    'ybinedges',yedges);
D = bwdist(~counts2);
D = -D;
figure; imagesc(D); axis square;
ids_im = watershed(D);
rgb = label2rgb(ids_im,'jet',[.5 .5 .5]);
figure; imagesc(rgb); axis square;
title('Watershed Transform')
watershed_ids = nan(N,1);
for ii = 1 : nbins
    ii_flags = ...
        pca_neuron_coeff(:,1) >= xedges(ii) & ...
        pca_neuron_coeff(:,1) < xedges(ii+1);
    for jj = 1 : nbins
        jj_flags = ...
            pca_neuron_coeff(:,2) >= yedges(jj) & ...
            pca_neuron_coeff(:,2) < yedges(jj+1);
        flags = ii_flags & jj_flags;
        watershed_ids(flags) = double(ids_im(ii,jj));
    end
end

%%
k_test = 4;
eval_mat = pca_neuron_coeff(:,1:2);
gm = fitgmdist(eval_mat,k_test);
% ids = kmeans(eval_mat,k_test);
% ids = cluster(gm,eval_mat);
ids = watershed_ids + 1; k_test = max(ids);
figure;
h = gscatter(eval_mat(:,1),eval_mat(:,2),ids);
A = nan(si_n_bins,k_test);
for kk = 1 : k_test
    A(:,kk) = nanmean(z(:,ids==kk),2);
end
Wa = pca(A);
[thetas,~] = cart2pol(Wa(:,1),Wa(:,2));
[~,theta_idcs] = sortrows(thetas);
theta_idcs = circshift(theta_idcs,sum(thetas<0));
figure; imagesc(A(:,theta_idcs)');
figure('position',[1.8000 41.8000 182.4000 740.8000],'color','w');
set(gca,'color','none','ycolor','none','nextplot','add','clipping','off');
for kk = 1 : k_test
    plot(si_time,A(:,theta_idcs(kk))+kk*3,...
        'color',h(theta_idcs(kk)).Color,...
        'linewidth',1.5);
    text(range(xlim)*-.1,kk*3,...
        sprintf('%i',sum(ids==theta_idcs(kk))));
end

%% cluster silhouettes
figure;
silhouette(pca_time_score(:,1:2),L);
silhouette(pca_neuron_coeff(:,1:2),ids);
% silhouette(pca_neuron_coeff(:,1:2),L2);
% silhouette(z',L2);
% silhouette(mon',L2);