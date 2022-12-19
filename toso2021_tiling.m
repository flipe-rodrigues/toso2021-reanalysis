%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct Ti-aligned, Ii-split psths
roi2use = [0,t_set(end)];
roi2plot = [-0,t_set(end)];
roi2use_n_bins = range(roi2use) * psthbin;
roi2plot_n_bins = range(roi2plot) * psthbin;
roi2use_time = linspace(roi2use(1),roi2use(2),roi2use_n_bins);
roi2plot_time = linspace(roi2plot(1),roi2plot(2),roi2plot_n_bins);
roi2use_flags = ...
    roi2plot_time >= roi2use(1) & ...
    roi2plot_time <= roi2use(2);
isi_n_bins = isi * psthbin;

% preallocation
s1_psths = nan(roi2plot_n_bins,n_neurons);
isi_psths = nan(isi_n_bins,n_neurons);
s2_psths = nan(roi2plot_n_bins,n_neurons);

% duration flags
t1_flags = t1 <= t_set(end);
t2_flags = t2 <= t_set(end);

% intensity flags
i1_flags = i1 == i_set(i1_mode_idx);
i2_flags = i2 == i_set(i2_mode_idx);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    s1_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i1_flags & ...
        t1_flags;
    isi_spike_flags = ...
        valid_flags & ...
        neuron_flags;
    s2_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i2_flags & ...
        t2_flags;
    if sum(s1_spike_flags) == 0 || ...
            sum(isi_spike_flags) == 0 || ...
            sum(s2_spike_flags) == 0
        continue;
    end
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%     s1_spike_flags = isi_spike_flags;
%     s2_spike_flags = isi_spike_flags;
    
    % fetch S1-aligned spike counts & compute spike rates
    s1_spike_counts = data.FR(s1_spike_flags,:);
    s1_spike_rates = conv2(...
        1,kernel.pdf,s1_spike_counts,'valid')' / psthbin * 1e3;
    s1_n_trials = size(s1_spike_counts,1);
    
    % fetch ISI-aligned spike counts & compute spike rates
    isi_spike_counts = data.FR(isi_spike_flags,:);
    isi_spike_rates = conv2(...
        1,kernel.pdf,isi_spike_counts,'valid')' / psthbin * 1e3;
    isi_n_trials = size(isi_spike_counts,1);
    
    % fetch S2-aligned spike counts & compute spike rates
    s2_spike_counts = data.FR(s2_spike_flags,:);
    s2_spike_rates = conv2(...
        1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
    s2_n_trials = size(s2_spike_counts,1);
    
    % S1-aligned spike rates
    s1_alignment_offset = ...
        pre_init_padding + ...
        pre_s1_delay(s1_spike_flags);
    s1_alignment_flags = ...
        valid_time >= s1_alignment_offset + roi2plot(1) & ...
        valid_time < s1_alignment_offset + t1(s1_spike_flags);
    s1_chunk_flags = ...
        valid_time >= s1_alignment_offset + roi2plot(1) & ...
        valid_time < s1_alignment_offset + roi2plot(2);
    s1_spkrates = s1_spike_rates;
    s1_spkrates(~s1_alignment_flags') = nan;
    s1_spkrates = reshape(...
        s1_spkrates(s1_chunk_flags'),...
        [roi2plot_n_bins,s1_n_trials])';
    
    % ISI-aligned spike rates
    isi_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(isi_spike_flags) + ...
        t1(isi_spike_flags);
    isi_alignment_flags = ...
        valid_time >= isi_alignment & ...
        valid_time < isi_alignment + isi;
    isi_chunk_flags = isi_alignment_flags;
    isi_spkrates = isi_spike_rates;
    isi_spkrates(~isi_alignment_flags') = nan;
    isi_spkrates = reshape(...
        isi_spkrates(isi_chunk_flags'),...
        [isi_n_bins,isi_n_trials])';
    
    % S2-aligned spike rates
    s2_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(s2_spike_flags) + ...
        t1(s2_spike_flags) + ...
        isi;
    s2_alignment_flags = ...
        valid_time >= s2_alignment + roi2plot(1) & ...
        valid_time < s2_alignment + t2(s2_spike_flags);
    s2_chunk_flags = ...
        valid_time >= s2_alignment + roi2plot(1) & ...
        valid_time < s2_alignment + roi2plot(2);
    s2_spkrates = s2_spike_rates;
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = reshape(...
        s2_spkrates(s2_chunk_flags'),...
        [roi2plot_n_bins,s2_n_trials])';
    
    % compute mean spike density functions
    s1_psths(:,nn) = nanmean(s1_spkrates,1);
    isi_psths(:,nn) = nanmean(isi_spkrates,1);
    s2_psths(:,nn) = nanmean(s2_spkrates,1);
end

%% normalization

% z-score S1-aligned spike density functions
s1_mus = nanmean(s1_psths(roi2use_flags,:),1);
s1_sigs = nanstd(s1_psths(roi2use_flags,:),0,1);
s1_zpsths = (s1_psths - s1_mus) ./ s1_sigs;

% z-score ISI-aligned spike density functions
isi_mus = nanmean(isi_psths,1);
isi_sigs = nanstd(isi_psths,0,1);
isi_zpsths = (isi_psths - isi_mus) ./ isi_sigs;

% z-score S2-aligned spike density functions
s2_mus = nanmean(s2_psths(roi2use_flags,:),1);
s2_sigs = nanstd(s2_psths(roi2use_flags,:),0,1);
s2_zpsths = (s2_psths - s2_mus) ./ s2_sigs;

%% PCA

% compute S1-aligned PCs
s1_coeff = pca(s1_zpsths(roi2use_flags,:));
s1_score = s1_zpsths * s1_coeff;

% compute ISI-aligned PCs
isi_coeff = pca(isi_zpsths);
isi_score = isi_zpsths * isi_coeff;

% compute S2-aligned PCs
weight_ref = median(squeeze(...
    surviving_trial_counts(flagged_neurons,i2_mode_idx,:)));

% compute observation weights
weights = ones(roi2use_n_bins,1) * max(weight_ref);
% for tt = 2 : n_t
%     t_flags = roi2use_time > t_set(tt-1) & roi2use_time <= t_set(tt);
%     weights(t_flags) = weight_ref(tt);
% end
weights = weights / max(weights);
s2_coeff = pca(s2_zpsths(roi2use_flags,:),...
    'weights',weights.^0);
% rica_mdl = rica(s2_zpsths,3);
% s2_coeff = rica_mdl.TransformWeights;
s2_score = s2_zpsths * s2_coeff;

%% colormap settings
clim = [-1.5,3.5];

%% S1-aligned tiling

% figure initialization
fig = figure(figopt,...
    ...'position',[1.8,41.8,516,740.8],...
    'name','tiling_s1');

% axes initialization
xxtick = unique([0;roi2use';roi2plot';t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
axes(axesopt.default,...
    ...'plotboxaspectratiomode','auto',...
    ...'ticklength',axesopt.default.ticklength * 2/3,...
    'xlim',roi2plot,...
    'ylim',[1,n_neurons],...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8));
title('S1-aligned PSTH raster');
xlabel('Time since S_1 onset (ms)');
ylabel('Neuron #');

% sort by angular position in PC space
[theta,~] = cart2pol(s1_coeff(:,1),s1_coeff(:,2));
[~,theta_idcs] = sortrows(theta);
theta_idcs = circshift(flipud(theta_idcs),125);

% plot selectivity heat map
imagesc(roi2plot,[1,n_neurons],s1_zpsths(:,theta_idcs)',clim);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = ...
    sprintf('z-score_{I_1 = %i %s}',...
    i_set(i1_mode_idx),i1_units);
clrlabel.string = 'Firing rate (z-score)';
clrlabel.fontsize = axesopt.default.fontsize * 1.1;
clrlabel.rotation = 270;
% clrlabel.position = [4.4,sum(clim)/2,0];
set(clrbar,...
    axesopt.colorbar,...
    'fontsize',axesopt.default.fontsize);
set(clrbar.Label,...
    clrlabel);

% plot alignment line
plot([1,1]*0,ylim,'--w');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% ISI-aligned tiling

% figure initialization
fig = figure(figopt,...
    ...'position',[1.8+516,41.8,516,740.8],...
    'name','tiling_isi');

% axes initialization
xxtick = 0:.5e3:isi;
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < isi) = {''};
axes(axesopt.default,...
    ...'plotboxaspectratiomode','auto',...
    ...'ticklength',axesopt.default.ticklength * 2/3,...
    'xlim',[0,isi],...
    'ylim',[1,n_neurons],...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8));
title('ISI-aligned PSTH raster');
xlabel('Time since S_1 offset / ISI onset (ms)');
ylabel('Neuron #');

% sort by angular position in PC space
[theta,~] = cart2pol(isi_coeff(:,1),isi_coeff(:,2));
[~,theta_idcs] = sortrows(theta);
theta_idcs = circshift(theta_idcs,-125);

% plot selectivity heat map
imagesc([0,isi],[1,n_neurons],isi_zpsths(:,theta_idcs)',clim);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = ...
    sprintf('z-score_{I_1 = %i %s}',...
    i_set(i1_mode_idx),i1_units);
clrlabel.string = 'Firing rate (z-score)';
clrlabel.fontsize = axesopt.default.fontsize * 1.1;
clrlabel.rotation = 270;
% clrlabel.position = [4.4,sum(clim)/2,0];
set(clrbar,...
    axesopt.colorbar,...
    'fontsize',axesopt.default.fontsize);
set(clrbar.Label,...
    clrlabel);

% plot alignment line
plot([1,1]*0,ylim,'--w');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% S2-aligned tiling

% figure initialization
fig = figure(figopt,...
    ...'position',[1.8+516*2,41.8,516,740.8],...
    'name','tiling_s2');

% axes initialization
xxtick = unique([0;roi2use';roi2plot';t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
axes(axesopt.default,...
    ...'plotboxaspectratiomode','auto',...
    ...'ticklength',axesopt.default.ticklength * 2/3,...
    'xlim',roi2plot,...
    'ylim',[1,n_neurons],...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8));
title('S2-aligned PSTH raster');
xlabel('Time since S_2 onset (ms)');
ylabel('Neuron #');

% sort by angular position in PC space
[theta,~] = cart2pol(s2_coeff(:,1),s2_coeff(:,2));
[~,theta_idcs] = sortrows(theta);
theta_idcs = circshift(theta_idcs,125);

% plot selectivity heat map
imagesc(roi2plot,[1,n_neurons],s2_zpsths(:,theta_idcs)',clim);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = sprintf('z-score_{I_2 = %i %s}',...
    i_set(i2_mode_idx),i2_units);
clrlabel.string = 'Firing rate (z-score)';
clrlabel.fontsize = axesopt.default.fontsize * 1.1;
clrlabel.rotation = 270;
% clrlabel.position = [4.4,sum(clim)/2,0];
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

%% S2-aligned PC projections
n_pcs2plot = 5;

% figure initialization
fig = figure(figopt,...
    'name','coeff_s2');

% axes initialization
sps = gobjects(n_pcs2plot,1);
for ii = 1 : n_pcs2plot
    sps(ii) = subplot(n_pcs2plot,1,ii);
    ylabel(sps(ii),sprintf('PC %i',ii));
end
xxtick = unique([0;roi2use';roi2plot';t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
set(sps,...
    axesopt.default,...
    'plotboxaspectratiomode','auto',...
    'xcolor','none',...
    'xlim',roi2plot,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel);
set(sps(end),...
    'xcolor','k');
xlabel(sps(end),'Time since S_2 onset (ms)');


% compute S2-aligned PCs
X = s2_psths(roi2use_flags,:);
% X(:,1:2:end) = ...
%     (2 * rand(1,ceil(n_neurons/2)) - 1) .* roi2use_time' + ...
%     randn(roi2use_n_bins,ceil(n_neurons/2)) * 10;
% X(:,2:2:end-1) = normpdf(roi2use_time,...
%     rand(floor(n_neurons/2),1) * 1e3, ...
%     25 + rand(floor(n_neurons/2),1) * 225)';
[Z,mu,sig] = zscore(X);
s2_coeff = pca(Z);
% ric = rica(Z,3);
% s2_coeff = ric.TransformWeights;
s2_score = Z * s2_coeff;
% sig = cov(X);
% [V,D] = eig(sig);
% V = rot90(V,2);
% D = rot90(D,2);
% S = X * V;

% iterate through PCs
for pc = 1 : n_pcs2plot

    % plot PC 1
    plot(sps(pc),roi2use_time,s2_score(:,pc),...
        'linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

% figure;
% hold on;
% ramp_flags = ismember(flagged_neurons,ramp_idcs.s1.on.up);
% plot(mean(s1_zpsths(:,ramp_flags),2))
% ramp_flags = ismember(flagged_neurons,ramp_idcs.s1.on.down);
% plot(mean(s1_zpsths(:,ramp_flags),2))
% 
% figure;
% hold on;
% ramp_flags = ismember(flagged_neurons,ramp_idcs.s2.on.up);
% plot(mean(s2_zpsths(:,ramp_flags),2))
% ramp_flags = ismember(flagged_neurons,ramp_idcs.s2.on.down);
% plot(mean(s2_zpsths(:,ramp_flags),2))

%% S2-aligned PC coefficient scatter

% figure initialization
fig = figure(figopt,...
    'name','coeff_s2');

% axes initialization
axes(axesopt.default,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'zlimspec','tight');
title('S2-aligned PC coefficients');
xlabel('PC 1');
ylabel('PC 2');

% compute S2-aligned PCs
X = s2_psths(roi2use_flags,:);
% X(:,1:2:end) = (2 * rand(1,ceil(n_neurons/2)) - 1) .* roi_time';
% X(:,2:2:end-1) = normpdf(roi_time,...
%     rand(floor(n_neurons/2),1) * 1e3, ...
%     rand(floor(n_neurons/2),1) * 250)';
% c = corr(X,repmat(roi_time,n_neurons,1)');
% [Z,mu,sig] = zscore(X);
% s2_coeff = pca(Z);
% sig = cov(X);
% [V,D] = eig(sig);
% V = rot90(V,2);
% D = rot90(D,2);
% S = X * V;

% plot selectivity heat map
% grapeplot(S(:,1),S(:,2));
% s2_coeff(:,1) = roi_time * Z;
% [~,s2_coeff(:,2)] = max(Z);
P = (Z - min(Z)) ./ range(Z);
P = P ./ nansum(P) + 1e-20;
h = -sum(P .* log2(P));
% s2_coeff(:,2) = h;
% s2_coeff(:,2) = var(diff(Z));
% 
% for nn = 1 : n_neurons
%     mdl = fitlm(Z(:,nn),roi_time);
%     s2_coeff(nn,1) = mdl.Coefficients.Estimate(2);
%     s2_coeff(nn,2) = mdl.Rsquared.Ordinary;
% end

% s2_coeff(:,1) = roi_time * Z;
% % P = (Z - min(Z)) ./ range(Z);
% s2_coeff(:,2) = h; % max(P); % var(diff(Z)); % sum(diff(P));

% ramp_flags = ismember(flagged_neurons,[ramp_idcs.s1.up;ramp_idcs.s1.down]);
s2_coeff = abs(s2_coeff);
n_hi = 5;
grapeplot(s2_coeff(:,1),sum(s2_coeff(:,2:n_hi),2));
% grapeplot(s2_coeff(ramp_flags,1),s2_coeff(ramp_flags,2),...
%     'markerfacecolor','b');
% grapeplot(s2_coeff(1:2:end,1),sum(s2_coeff(1:2:end,2:n_hi),2),...
%     'markerfacecolor','b');

% update axes limits
xxlim = xlim;
yylim = ylim;
zzlim = zlim;
set(gca,...
    'xlim',xxlim+[-1,1]*.05*range(xxlim),...
    'ylim',yylim+[-1,1]*.05*range(yylim),...
    'zlim',zzlim+[-1,1]*.05*range(zzlim),...
    'xtick',unique([0,xxlim]),...
    'ytick',unique([0,yylim]),...
    'ztick',unique([0,zzlim]),...
    'xticklabel',{'','0',''},...
    'yticklabel',{'','0',''},...
    'zticklabel',{'','0',''});

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end