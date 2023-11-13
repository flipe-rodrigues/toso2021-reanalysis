%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct T1-offset-aligned, contrast-split psths
pre_padd = 500;
roi2use = [0,isi];
roi2plot = [-t_set(end),isi+t_set(end)];
roi2use_n_bins = range(roi2use) * psthbin;
roi2plot_n_bins = range(roi2plot) * psthbin;
roi2use_time = linspace(roi2use(1),roi2use(2),roi2use_n_bins);
roi2plot_time = linspace(roi2plot(1),roi2plot(2),roi2plot_n_bins);
roi2use_flags = ...
    roi2plot_time >= roi2use(1) & ...
    roi2plot_time <= roi2use(2);

% preallocation
psths = nan(roi2plot_n_bins,n_neurons,n_contrasts);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        spike_flags = ...
            valid_flags & ...
            neuron_flags & ...
            contrast_flags;
        if sum(spike_flags) == 0
            continue;
        end
        
        % fetch T1-offset-aligned spike counts & compute spike rates
        spike_counts = data.FR(spike_flags,:);
        spike_rates = conv2(...
            1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
        n_trials = size(spike_counts,1);
        
        % T1-offset-aligned spike rates
        alignment_offset = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags);
        alignment_flags = ...
            valid_time >= alignment_offset - t1(spike_flags) & ...
            valid_time < alignment_offset + isi + t2(spike_flags);
        chunk_flags = ...
            valid_time >= alignment_offset + roi2plot(1) & ...
            valid_time < alignment_offset + roi2plot(2);
        spkrates = spike_rates;
        spkrates(~alignment_flags') = nan;
        spkrates = reshape(...
            spkrates(chunk_flags'),...
            [roi2plot_n_bins,n_trials])';
        
        % compute mean spike density functions
        psths(:,nn,ii) = nanmean(spkrates,1);
    end
end

% nan handling
% psths(isnan(psths)) = 0;

%% normalization

% z-score T1-offset-aligned spike density functions
mus = nanmean(psths,[1,3]);
sigs = nanstd(psths,0,[1,3]);
zpsths = (psths - mus) ./ sigs;

%% cross-condition concatenations

% concatenate T1-offset-aligned psths across conditions
concat_all = nan(roi2use_n_bins*n_contrasts,n_neurons);
concat_extr = nan(roi2use_n_bins*(n_contrasts-1),n_neurons);
concat_mode = nan(roi2use_n_bins,n_neurons);
concat_diff = nan(roi2use_n_bins,n_neurons);
for nn = 1 : n_neurons
    nn_zpsths_all = zpsths(roi2use_flags,nn,:);
    nn_zpsths_extr = zpsths(roi2use_flags,nn,(1:n_contrasts)~=contrast_mode_idx);
    nn_zpsths_mode = zpsths(roi2use_flags,nn,contrast_mode_idx);
    nn_zpsths_diff = zpsths(roi2use_flags,nn,end) - zpsths(roi2use_flags,nn,1);
    concat_all(:,nn) = nn_zpsths_all(:);
    concat_extr(:,nn) = nn_zpsths_extr(:);
    concat_mode(:,nn) = nn_zpsths_mode(:);
    concat_diff(:,nn) = nn_zpsths_diff(:);
end

%% PCA
[coeff,~,~,~,exp] = pca(concat_mode);
exp_cutoff = 80;
n_pcs2use = max(sum(cumsum(exp) <= exp_cutoff),2);

%% compute entropy
fudge = 1e-6;
p = concat_mode - min(concat_mode);
p = p ./ nansum(p) + fudge;
ent = -sum(p.*log2(p));
[~,ent_idcs] = sort(ent);

%% angular position in pc space
[theta,rho] = cart2pol(coeff(:,1),coeff(:,2));
[~,theta_idcs] = sortrows(theta);
% theta_idcs = circshift(theta_idcs,sum(theta>0));
theta_idcs = circshift(theta_idcs,-80);

%% agglomerative hierarchical clustering
diss = pdist(concat_mode','euclidean');
diss = pdist(coeff(:,1:n_pcs2use),'euclidean');
% diss = pdist([coeff(:,1:n_pcs2use),ent',theta],'euclidean');
tree = linkage(diss,'ward');
leaforder = optimalleaforder(tree,diss);

%% dendrogram

% figure initialization
figure(figopt,...
    'position',[1.0554e+03 905 359.2000 1.3176e+03]);

% plot dendrogram
clr_threshold = max(tree(:,3)) * 2/3;
[d,T,leaf_idcs] = dendrogram(tree,0,...
    'orientation','right',...
    'reorder',leaforder,...
    'colorthreshold',clr_threshold);

% update axes
set(gca,...
    axesopt.default,...
    'ylim',[1,n_neurons],...
    'plotboxaspectratiomode','auto');
xlabel('Euclidean Distance');
ylabel('Neuron #');

% plot color threshold
plot([1,1]*clr_threshold,ylim,'--k');

%% PC weights

% figure initialization
figure(figopt,...
    'position',[903.4000 905 150 1.3176e+03]);
axes(axesopt.default,...
    'xlim',[1,n_pcs2use]+[-1,1]*.5,...
    'ylim',[1,n_neurons],...
    'xtick',1:n_pcs2use,...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8),...
    'plotboxaspectratiomode','auto');
xlabel('PC #');
ylabel('Neuron #');

% plot selectivity heat map
imagesc([1,n_pcs2use],[1,n_neurons],coeff(leaf_idcs,1:n_pcs2use));

%% heatmap

% figure initialization
figure(figopt,...
    'position',[541.8000 905 359.2000 1.3176e+03]);
axes(axesopt.default,...
    'xlim',roi2plot,...
    'ylim',[1,n_neurons],...
    'xtick',unique([0,roi2use,roi2plot]),...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8),...
    'plotboxaspectratiomode','auto');
xlabel('Time since T_1 offset (ms)');
ylabel('Neuron #');

% plot selectivity heat map
clim = [-1,1] * 3;
imagesc(roi2plot,[1,n_neurons],...
    squeeze(zpsths(:,leaf_idcs,contrast_mode_idx))',clim);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = ...
    sprintf('z-score_{%s = %i %s}',...
    contrast_lbl,contrast_set(contrast_mode_idx),contrast_units);
clrlabel.fontsize = axesopt.default.fontsize * 1.1;
clrlabel.rotation = 270;
clrlabel.position = [4.4,0,0];
clrlabel.position = [4.4,sum(clim)/2,0];
set(clrbar,...
    axesopt.colorbar,...
    'fontsize',axesopt.default.fontsize);
set(clrbar.Label,...
    clrlabel);

% plot alignment line
plot([1,1]*0,ylim,'--w');

%% average psths per cluster
max_n_clusters = 6;
for cc = 1 : max_n_clusters
    
    % figure initialization
    figure(figopt,...
        'windowstyle','docked',...
        'name',sprintf('cluster averages: %i clusters',cc));
    
    % axes initialization
    sps = gobjects(cc,1);
    for kk = 1 : cc
        sps(kk) = subplot(cc,1,kk);
        ylabel(sps(kk),'Firing rate (z-score)');
    end
    xlabel(sps(end),'Time since T_1 offset (ms)');
    set(sps,...
        axesopt.default,...
        'xlim',roi2plot,...
        'xtick',unique([0,roi2use,roi2plot]),...
        'plotboxaspectratiomode','auto',...
        'layer','top');
    
    % clustering
    T = cluster(tree,...
        'maxclust',cc,...
        'criterion','distance');
    
    % iterate through clusters
    for kk = 1 : cc
        cluster_flags = T == kk;
        
        % iterate through contrasts
        for ii = 1 : n_contrasts
            
            % compute cluster stats
            mu = nanmean(zpsths(:,cluster_flags,ii),2)';
            sig = nanstd(zpsths(:,cluster_flags,ii),0,2)';
            sem = sig / sqrt(sum(cluster_flags));
            
            % patch s.e.m.
            nan_flags = isnan(mu);
            xpatch = [roi2plot_time(~nan_flags),fliplr(roi2plot_time(~nan_flags))];
            ypatch = [mu(~nan_flags)-sem(~nan_flags),fliplr(mu(~nan_flags)+sem(~nan_flags))];
            patch(sps(kk),xpatch,ypatch,contrast_clrs(ii,:),...
                'facealpha',.25,...
                'edgecolor','none');
            
            % plot average through time
            plot(sps(kk),roi2plot_time,mu,...
                'linewidth',1.5,...
                'color',contrast_clrs(ii,:));
            
            % plot average at T1-onset
            onset_flags = [true,nan_flags(1:end-1)] & ~nan_flags;
            plot(sps(kk),roi2plot_time(onset_flags),mu(onset_flags),...
                'linewidth',1.5,...
                'marker','o',...
                'markersize',7.5,...
                'markerfacecolor','w',...
                'markeredgecolor',contrast_clrs(ii,:));
            
            % highlight contrast mode
            if ii == contrast_mode_idx
                
                % plot average through time
                roi2use_flags = ...
                    roi2plot_time >= roi2use(1) & ...
                    roi2plot_time <= roi2use(2);
                p = plot(sps(kk),roi2use_time,mu(roi2use_flags),...
                    'linewidth',3,...
                    'linestyle','-',...
                    'color','k');
            end
        end
        
        % plot alignment line
        yylim = ylim(sps(kk));
        plot(sps(kk),[1,1]*0,yylim,'--k');
        plot(sps(kk),[1,1]*isi,yylim,'--k');
        ylim(sps(kk),yylim);
        
        % text annotation
        text(sps(kk),.05,.95,sprintf('N = %i neurons',sum(cluster_flags)),...
            'units','normalized');
        
        % ui restacking
        uistack(p,'top');
    end
end