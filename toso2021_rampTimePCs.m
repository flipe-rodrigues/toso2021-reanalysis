%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% Gap statistic settings
max_k = 6;
max_pc2consider = 10;
pc_clrs = flipud(gray(max_pc2consider+1));
pc_clrs = pc_clrs(2:end,:);

%% ROI settings
si_roi = [0,t_set(end)];
si_n_bins = range(si_roi) / psthbin;
si_time = linspace(si_roi(1),si_roi(2),si_n_bins);

%% construct Si- and ISI-aligned psths

% preallocation
s1_psths = nan(si_n_bins,n_neurons);
s2_psths = nan(si_n_bins,n_neurons);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    spike_flags = ...
        valid_flags & ...
        neuron_flags;
    
    % fetch spike rates
    spike_rates = data.SDF(spike_flags,:);
    n_trials = size(spike_rates,1);
    
    % S1-aligned spike rates
    s1_alignment_offset = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags);
    s1_alignment_flags = ...
        padded_time >= s1_alignment_offset + si_roi(1) & ...
        padded_time < s1_alignment_offset + t1(spike_flags);
    s1_chunk_flags = ...
        padded_time >= s1_alignment_offset + si_roi(1) & ...
        padded_time < s1_alignment_offset + si_roi(2);
    s1_spkrates = spike_rates';
    s1_spkrates(~s1_alignment_flags') = nan;
    s1_spkrates = reshape(...
        s1_spkrates(s1_chunk_flags'),...
        [si_n_bins,n_trials])';
    
    % S2-aligned spike rates
    s2_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(spike_flags) + ...
        t1(spike_flags) + ...
        isi;
    s2_alignment_flags = ...
        padded_time >= s2_alignment + si_roi(1) & ...
        padded_time < s2_alignment + t2(spike_flags);
    s2_chunk_flags = ...
        padded_time >= s2_alignment + si_roi(1) & ...
        padded_time < s2_alignment + si_roi(2);
    s2_spkrates = spike_rates';
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = reshape(...
        s2_spkrates(s2_chunk_flags'),...
        [si_n_bins,n_trials])';
    
    % compute mean spike density functions
    s1_psths(:,nn) = nanmean(s1_spkrates);
    s2_psths(:,nn) = nanmean(s2_spkrates);
end

%% normalization

% z-score S1-aligned spike density functions
s1_mus = nanmean(s1_psths,[1,3]);
s1_sigs = nanstd(s1_psths,0,[1,3]);
s1_zpsths = (s1_psths - s1_mus) ./ s1_sigs;

% z-score S2-aligned spike density functions
s2_mus = nanmean(s2_psths,[1,3]);
s2_sigs = nanstd(s2_psths,0,[1,3]);
s2_zpsths = (s2_psths - s2_mus) ./ s2_sigs;

%% compute Si-aligned weights
time_mat = repmat(si_time,n_total_trials,1);
s1_weights = sum(time_mat(valid_flags,:) <= t1(valid_flags));
s1_weights = s1_weights / sum(s1_weights);
s2_weights = sum(time_mat(valid_flags,:) <= t2(valid_flags));
s2_weights = s2_weights / sum(s2_weights);

%% PCA

% compute S1-aligned PCs
[s1_coeff,s1_score,~,~,s1_exp] = pca(zscore(s1_psths'));

% compute S2-aligned PCs
[s2_coeff,s2_score,~,~,s2_exp] = pca(zscore(s2_psths'));

%% compute gap statistics

% preallocation
s1_eval = struct(...
    'gap',nan(max_pc2consider,max_k),...
    'se',nan(max_pc2consider,max_k),...
    'k',nan(max_pc2consider,1));
s2_eval = struct(...
    'gap',nan(max_pc2consider,max_k),...
    'se',nan(max_pc2consider,max_k),...
    'k',nan(max_pc2consider,1));

% iterate through number of PCs considered
for pc = 1 : max_pc2consider
    progressreport(pc,max_pc2consider,'computing gap statistics');
    
    % evaluate S1 coefficients
    eval = evalclusters(s1_coeff(:,1:pc),'kmeans','gap',...
        'referencedistribution','uniform',...
        'klist',1:max_k);
    s1_eval.gap(pc,:) = eval.CriterionValues;
    s1_eval.se(pc,:) = eval.SE;
    s1_eval.k(pc) = eval.OptimalK;
    
    % evaluate S2 coefficients
    eval = evalclusters(s2_coeff(:,1:pc),'kmeans','gap',...
        'referencedistribution','uniform',...
        'klist',1:max_k);
    s2_eval.gap(pc,:) = eval.CriterionValues;
    s2_eval.se(pc,:) = eval.SE;
    s2_eval.k(pc) = eval.OptimalK;
end

%% parse cluster indices & colors

% S1-aligned clusters
s1_ramp_flags = ismember(flagged_neurons,cluster_idcs.s1{'ramp'});
s1_nonramp_flags = ismember(flagged_neurons,cluster_idcs.s1{'nonramp'});
s1_cluster_ids = s1_ramp_flags + s1_nonramp_flags * 2;
s1_cluster_flags = s1_cluster_ids > 0;
s1_cluster_clrs = repmat([1,1,1],n_neurons,1);
s1_cluster_clrs(s1_cluster_flags,:) = ...
    ramp_clrs(s1_cluster_ids(s1_cluster_flags),:);

% S2-aligned clusters
s2_ramp_flags = ismember(flagged_neurons,cluster_idcs.s2{'ramp'});
s2_nonramp_flags = ismember(flagged_neurons,cluster_idcs.s2{'nonramp'});
s2_cluster_ids = s2_ramp_flags + s2_nonramp_flags * 2;
s2_cluster_flags = s2_cluster_ids > 0;
s2_cluster_clrs = repmat([1,1,1],n_neurons,1);
s2_cluster_clrs(s2_cluster_flags,:) = ...
    ramp_clrs(s2_cluster_ids(s2_cluster_flags),:);

%% S1-aligned tiling

% figure initialization
fig = figure(figopt,...
    'position',[100,200,560,420],...
    'name','tiling_s1');

% axes initialization
xxtick = unique([0;si_roi';t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
axes(axesopt.default,...
    'xlim',si_roi,...
    'ylim',[1,n_neurons],...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8),...
    'clipping','off');
title('S1-aligned PSTH raster');
xlabel('Time since S_1 onset (ms)');
ylabel({'Neuron #','(sorted by S1-aligned PCs)'});

% sort by angular position in PC space
[s1_theta,~] = cart2pol(s1_coeff(:,1),s1_coeff(:,2));
[~,s1_theta_idcs] = sortrows(s1_theta);
s1_theta_idcs = circshift(s1_theta_idcs,150);

% color limits
clim = [-2,4];

% plot psth raster
imagesc(si_roi,[1,n_neurons],s1_zpsths(:,s1_theta_idcs)',clim);

% plot cluster affordances
plot(0-.025*range(xlim),find(s1_ramp_flags(s1_theta_idcs)),...
    'color',ramp_clrs(1,:),...
    'marker','.',...
    'markersize',5);
plot(0-.0125*range(xlim),find(s1_nonramp_flags(s1_theta_idcs)),...
    'color',ramp_clrs(2,:),...
    'marker','.',...
    'markersize',5);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = sprintf('z-score_{I_2 = %i %s}',...
    i_set(i2_mode_idx),i2_units);
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

% plot alignment line
plot([1,1]*0,ylim,'--w',...
    'linewidth',1.5);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% S2-aligned tiling

% figure initialization
fig = figure(figopt,...
    'position',[100,200,560,420],...
    'name','tiling_s2');

% axes initialization
xxtick = unique([0;si_roi';t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
axes(axesopt.default,...
    'xlim',si_roi,...
    'ylim',[1,n_neurons],...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8),...
    'clipping','off');
title('S2-aligned PSTH raster');
xlabel('Time since S_2 onset (ms)');
ylabel({'Neuron #','(sorted by S2-aligned PCs)'});

% sort by angular position in PC space
[s2_theta,~] = cart2pol(s2_coeff(:,1),s2_coeff(:,2));
[~,s2_theta_idcs] = sortrows(s2_theta);
s2_theta_idcs = flipud(circshift(s2_theta_idcs,-125));

% plot psth raster
imagesc(si_roi,[1,n_neurons],s2_zpsths(:,s2_theta_idcs)',clim);

% plot cluster affordances
plot(0-.025*range(xlim),find(s2_ramp_flags(s2_theta_idcs)),...
    'color',ramp_clrs(1,:),...
    'marker','.',...
    'markersize',5);
plot(0-.0125*range(xlim),find(s2_nonramp_flags(s2_theta_idcs)),...
    'color',ramp_clrs(2,:),...
    'marker','.',...
    'markersize',5);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = sprintf('z-score_{I_2 = %i %s}',...
    i_set(i2_mode_idx),i2_units);
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

% plot alignment line
plot([1,1]*0,ylim,'--w',...
    'linewidth',1.5);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% S1-aligned PC coefficient scatter
fig = figure(figopt,...
    'name','pc_coefficients2D_s1');
axes(axesopt.default,...
    'clipping','off');
xlabel('S_{1}-aligned PC 1 coefficient_{1}');
ylabel('S_{1}-aligned PC 2 coefficient_{2}');

% coefficient scatter
grapeplot(s1_score(s2_cluster_flags,1),s1_score(s2_cluster_flags,2),...
    'markerfacecolor',s1_cluster_clrs(s2_cluster_flags,:));

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

%% S2-aligned PC coefficient scatter
fig = figure(figopt,...
    'name','pc_coefficients2D_s2');
axes(axesopt.default,...
    'clipping','off');
xlabel('S_{2}-aligned PC 1 coefficient_{1}');
ylabel('S_{2}-aligned PC 2 coefficient_{2}');

% coefficient scatter
grapeplot(s2_score(s2_cluster_flags,1),s2_score(s2_cluster_flags,2),...
    'markerfacecolor',s2_cluster_clrs(s2_cluster_flags,:));

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

%% plot S1-aligned gap statistics
fig = figure(figopt,...
    'name','pc_gap_s1');
axes(axesopt.default,...
    'xlim',[1,max_k]+[-1,1]*.05*max_k,...
    'xtick',1:max_k,...
    'ylimspec','tight');
xlabel('Number of clusters');
ylabel('Gap value');

% iterate through number of PCs considered
for pc = 1 : max_pc2consider
    errorbar(1:max_k,s1_eval.gap(pc,:),s1_eval.se(pc,:),...
        'color',pc_clrs(pc,:),...
        'linewidth',1.5,...
        'capsize',0);
    plot(s1_eval.k(pc),s1_eval.gap(pc,s1_eval.k(pc)),...
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
    @(x)sprintf('first %i PCs_{s1}',x),1:max_pc2consider,...
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

%% plot S2-aligned gap statistics
fig = figure(figopt,...
    'name','pc_gap_s2');
axes(axesopt.default,...
    'xlim',[1,max_k]+[-1,1]*.05*max_k,...
    'xtick',1:max_k,...
    'ylimspec','tight');
xlabel('Number of clusters');
ylabel('Gap value');

% iterate through number of PCs considered
for pc = 1 : max_pc2consider
    errorbar(1:max_k,s2_eval.gap(pc,:),s2_eval.se(pc,:),...
        'color',pc_clrs(pc,:),...
        'linewidth',1.5,...
        'capsize',0);
    plot(s2_eval.k(pc),s2_eval.gap(pc,s2_eval.k(pc)),...
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
    @(x)sprintf('first %i PCs_{s2}',x),1:max_pc2consider,...
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

%%
load('r.mat','r_bumps');
r_ramps = si_time' .* randn(1,size(r,2)) + randn(si_n_bins,size(r,2))*0;
figure; plot(si_time,r_ramps(:,1:10),'k');
figure; imagesc(r_ramps');
r_comb = [r_ramps,r];
z_comb = zscore(r_comb);
figure; imagesc(z_comb');

%%
max_k = 15;
fig = figure(figopt,...
    'name','pc_coefficients_gap_s2');
axes(axesopt.default,...
    'xlim',[1,max_k]+[-1,1]*.05*max_k,...
    'xtick',1:max_k,...
    'ylimspec','tight');
xlabel('Number of clusters');
ylabel('Gap value');

Z = zscore(s2_psths');
% Z = z_comb;
[~,W] = pca(Z);
% t = si_time;
% for ii = 1 : size(Z,2)
%     W(ii,1) = Z(:,ii) \ t';
%     p = polyfit(t,Z(:,ii),2);
% %     W(ii,1) = p(3);
%     z = Z(:,ii);
%     dz = diff(z);
%     W(ii,2) = entropy(diff(Z(:,ii)));
%     W(ii,2) = p(3);
%     %      z = Z(:,ii) - min(Z(:,ii));
%     %      W(ii,2) = t * (z / nansum(z));
% end

% whiten
% W = W';
% sig = (W' * W) / n_neurons;
% [U,S] = svd(sig);
% tol = eps(class(W));
% W = W - mean(W,1);
% W = W * U * diag(1 ./ sqrt(diag(S) + tol)) * U';
% W = W';

eval = evalclusters(W(:,1:2),'kmeans','gap',...
    'referencedistribution','uniform',...
    'klist',1:max_k);
s2_eval2.gap = eval.CriterionValues;
s2_eval2.se = eval.SE;
s2_eval2.k = eval.OptimalK;

% iterate through number of PCs considered
errorbar(1:max_k,s2_eval2.gap,s2_eval2.se,...
    'color','k',...
    'linewidth',1.5,...
    'capsize',0);
plot(s2_eval2.k,s2_eval2.gap(s2_eval2.k),...
        'color','k',...
        'marker','o',...
        'markersize',7.5,...
        'markeredgecolor','k',...
        'markerfacecolor','w',...
        'linewidth',1.5,...
        'handlevisibility','off');
    
% update axis
set(gca,...
    'ytick',ylim,...
    'ylim',ylim+[-1,1]*.05*range(ylim),...
    'yticklabel',{});

%%
figure; h = gscatter(W(:,1),W(:,2),eval.OptimalY);
A = nan(si_n_bins,eval.OptimalK);
for kk = 1 : eval.OptimalK
    A(:,kk) = nanmean(Z(:,eval.OptimalY==kk),2);
end
Wa = pca(A);
[theta,~] = cart2pol(Wa(:,1),Wa(:,2));
[~,theta_idcs] = sortrows(theta);
theta_idcs = circshift(theta_idcs,sum(theta<0));
figure; imagesc(A(:,theta_idcs)');
figure('position',[1.8000 41.8000 182.4000 740.8000],'color','w');
set(gca,'color','none','ycolor','none','nextplot','add','clipping','off');
for kk = 1 : eval.OptimalK
    plot(si_time,A(:,theta_idcs(kk))+kk*3,...
        'color',h(theta_idcs(kk)).Color,...
        'linewidth',1.5);
    text(range(xlim)*-.1,kk*3,...
        sprintf('%i',sum(eval.OptimalY==theta_idcs(kk))));
end