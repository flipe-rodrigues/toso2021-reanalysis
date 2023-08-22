%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct Ti-aligned, Ii-split psths
pre_padd = 500;
roi2use = [-0,t_set(end)];
roi2plot = [-pre_padd,t_set(end)];
roi2plot_padded = roi2plot + [-1,1] * .05 * range(roi2plot);
roi2use_n_bins = range(roi2use) / psthbin;
roi2plot_n_bins = range(roi2plot_padded) / psthbin;
roi2use_time = linspace(roi2use(1),roi2use(2),roi2use_n_bins);
roi2plot_time = linspace(roi2plot_padded(1),roi2plot_padded(2),roi2plot_n_bins);
roi2use_flags = ...
    roi2plot_time >= roi2use(1) & ...
    roi2plot_time <= roi2use(2);
isi_n_bins = isi / psthbin;

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
        ...i1_flags & ...
        t1_flags;
    isi_spike_flags = ...
        valid_flags & ...
        neuron_flags;
    s2_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        ...i2_flags & ...
        t2_flags;
    if sum(s1_spike_flags) == 0 || ...
            sum(isi_spike_flags) == 0 || ...
            sum(s2_spike_flags) == 0
        continue;
    end

    % fetch S1-aligned spike counts & compute spike rates
    s1_spike_counts = data.FR(s1_spike_flags,:);
    s1_spike_rates = data.SDF(s1_spike_flags,:);
    s1_n_trials = size(s1_spike_counts,1);
    
    % fetch ISI-aligned spike counts & compute spike rates
    isi_spike_counts = data.FR(isi_spike_flags,:);
    isi_spike_rates = data.SDF(isi_spike_flags,:);
    isi_n_trials = size(isi_spike_counts,1);
    
    % fetch S2-aligned spike counts & compute spike rates
    s2_spike_counts = data.FR(s2_spike_flags,:);
    s2_spike_rates = data.SDF(s2_spike_flags,:);
    s2_n_trials = size(s2_spike_counts,1);
    
    % S1-aligned spike rates
    s1_alignment_offset = ...
        pre_init_padding + ...
        pre_s1_delay(s1_spike_flags);
    s1_alignment_flags = ...
        valid_time >= s1_alignment_offset + roi2plot_padded(1) & ...
        valid_time < s1_alignment_offset + t1(s1_spike_flags);
    s1_chunk_flags = ...
        valid_time >= s1_alignment_offset + roi2plot_padded(1) & ...
        valid_time < s1_alignment_offset + roi2plot_padded(2);
    s1_spkrates = s1_spike_rates';
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
    isi_spkrates = isi_spike_rates';
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
        valid_time >= s2_alignment + roi2plot_padded(1) & ...
        valid_time < s2_alignment + t2(s2_spike_flags);
    s2_chunk_flags = ...
        valid_time >= s2_alignment + roi2plot_padded(1) & ...
        valid_time < s2_alignment + roi2plot_padded(2);
    s2_spkrates = s2_spike_rates';
    s2_spkrates(~s2_alignment_flags') = nan;
    s2_spkrates = reshape(...
        s2_spkrates(s2_chunk_flags'),...
        [roi2plot_n_bins,s2_n_trials])';
    
    % compute mean spike density functions
    s1_psths(:,nn) = nanmean(s1_spkrates,1);
    isi_psths(:,nn) = nanmean(isi_spkrates,1);
    s2_psths(:,nn) = nanmean(s2_spkrates,1);
end

%% compute alignment-specific weights

% compute S1-aligned weights
time_mat = repmat(roi2use(1) + psthbin : psthbin : roi2use(2),n_total_trials,1);
s1_weights = sum(time_mat(valid_flags,:) <= t1(valid_flags));
s1_weights = s1_weights / sum(s1_weights);

% compute ISI-aligned weights
isi_weights = ones(1,isi_n_bins) / isi_n_bins;

% compute S2-aligned weights
time_mat = repmat(roi2use(1) + psthbin : psthbin : roi2use(2),n_total_trials,1);
s2_weights = sum(time_mat(valid_flags,:) <= t2(valid_flags));
s2_weights = s2_weights / sum(s2_weights);

%% normalization

% z-score S1-aligned spike density functions
s1_mus = nanmean(s1_psths(roi2use_flags,:),1);
s1_sigs = nanstd(s1_psths(roi2use_flags,:),0,1);
% s1_mus = s1_weights(roi2use_flags) * s1_psths(roi2use_flags,:);
% s1_sigs = nanstd(s1_psths(roi2use_flags,:),s1_weights(roi2use_flags),1);
s1_zpsths = (s1_psths - s1_mus) ./ s1_sigs;

% z-score ISI-aligned spike density functions
isi_mus = nanmean(isi_psths,1);
isi_sigs = nanstd(isi_psths,0,1);
isi_zpsths = (isi_psths - isi_mus) ./ isi_sigs;

% z-score S2-aligned spike density functions
s2_mus = nanmean(s2_psths(roi2use_flags,:),1);
s2_sigs = nanstd(s2_psths(roi2use_flags,:),0,1);
% s2_mus = s2_weights(roi2use_flags) * s2_psths(roi2use_flags,:);
% s2_sigs = nanstd(s2_psths(roi2use_flags,:),s2_weights(roi2use_flags),1);
s2_zpsths = (s2_psths - s2_mus) ./ s2_sigs;

%% PCA

% compute S1-aligned PCs
s1_coeff = pca(s1_zpsths(roi2use_flags,:),...
    'weights',s1_weights);
s1_score = s1_zpsths * s1_coeff;

% compute ISI-aligned PCs
isi_coeff = pca(isi_zpsths,...
    'weights',isi_weights);
isi_score = isi_zpsths * isi_coeff;

% compute S2-aligned PCs
s2_coeff = pca(s2_zpsths(roi2use_flags,:),...
    'weights',s2_weights);
s2_score = s2_zpsths * s2_coeff;

%%

for pc_s1 = 1 : 3
    figure('position',[94.6000+(pc_s1-1)*560 41.8000 560 1.0288e+03]);
    [x,idcs] = sort(s1_coeff(:,pc_s1));
    for pc_s2 = 1 : 3
        subplot(3,1,pc_s2);
        hold on;
        xlabel(sprintf('S1-aligned PC_%i weights',pc_s1))
        ylabel(sprintf('S2-aligned PC_%i weights',pc_s2))
        y = s2_coeff(idcs,pc_s2);
        mdl = fitlm(x,y);
        if mdl.Coefficients.pValue(end) < .05
            linestyle = '-';
        else 
            linestyle = '--';
        end
        plot(x,y,'.k')
        plot(x,mdl.predict(x),'k',...
            'linestyle',linestyle,...
            'linewidth',1.5);
        axis tight;
        axis square;
    end
end

%%

figure;
hold on;
xlabel('S1-aligned variance')
ylabel('S2-aligned variance')
[x,idcs] = sort(nanvar(s1_psths));
y = nanvar(s2_psths(:,idcs));
mdl = fitlm(x,y);
if mdl.Coefficients.pValue(end) < .05
    linestyle = '-';
else
    linestyle = '--';
end
plot(x,y,'.k')
plot(x,mdl.predict(x'),'k',...
    'linestyle',linestyle,...
    'linewidth',1.5);
axis tight;
axis square;

plot(xlim,xlim,':')

set(gca,...
    'xlim',quantile(x,[0,.95]),...
    'ylim',quantile(y,[0,.95]),...
    'xscale','linear',...
    'yscale','linear');

%% colormap settings
clim = [-1.5,3.5];

%% S1-aligned tiling

% figure initialization
fig = figure(figopt,...
    'position',[10 630 560 420],...[1.8,41.8,516,740.8],...
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
theta_idcs = circshift(theta_idcs,-200);

% plot psth raster
imagesc(roi2plot_padded,[1,n_neurons],s1_zpsths(:,theta_idcs)',clim);

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
    'position',[440 630 560*2 420],...[1.8+516,41.8,516,740.8],...
    'name','tiling_isi');

% axes initialization
xxtick = 0:.5e3:isi;
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < isi) = {''};
axes(axesopt.default,...
    'plotboxaspectratio',[2,1,1],...
    'ticklength',axesopt.default.ticklength * 1/2,...
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
% theta_idcs = circshift(theta_idcs,-125);

% plot psth raster
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
    'position',[1325 630 560 420],...[1.8+516*2,41.8,516,740.8],...
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
theta_idcs = circshift(theta_idcs,40);

% plot psth raster
imagesc(roi2plot_padded,[1,n_neurons],s2_zpsths(:,theta_idcs)',clim);

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
return;

%% PC projections

% figure initialization
fig = figure(figopt,...
    'position',[535,130,966,860],...
    'name',sprintf('pc_projections_s2_%s',contrast_str));
n_pcs2plot = 6;
sps = gobjects(n_pcs2plot,1);
for pc = 1 : n_pcs2plot
    sp_idx = pc * 2 - 1 - (pc > n_pcs2plot / 2) * (n_pcs2plot - 1);
    sps(pc) = subplot(n_pcs2plot/2,2,sp_idx);
    xlabel(sps(pc),'Time since S_2 onset (ms)');
    %     ylabel(sps(pc),sprintf('PC %i\n%.1f%% variance',pc,exp_pca(pc)));
    ylabel(sps(pc),sprintf('PC %i',pc));
end
xxtick = unique([roi2plot';0;t_set]);
xxticklabel = num2cell(xxtick);
xxticklabel(xxtick > 0 & xxtick < t_set(end)) = {''};
set(sps,...
    axesopt.default,...
    'xlim',roi2plot + [-1,1] * .05 * range(roi2plot),...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'plotboxaspectratio',[3,1,1]);

% link axes
linkaxes(sps,'x');

% iterate through pcs
for pc = 1 : n_pcs2plot
    trial_flags = ...
        valid_flags;
    
    % compute surviving trial counts (through time)
    time_mat = repmat(...
        padded_roi(1) + psthbin : psthbin : padded_roi(2),n_total_trials,1);
    surviving_trial_counts = sum(time_mat(trial_flags,:) <= t2(trial_flags));
    
    % patch projection
    mu_xpatch = roi2plot_time;
    mu_ypatch = s2_score(:,pc)';
    mu_ypatch(end) = nan;
    mu_apatch = surviving_trial_counts;
    mu_apatch = mu_apatch ./ max(mu_apatch) .* ...
        range(alphabounds_mu) + alphabounds_mu(1);
    if fadeifnoisy
        alpha_levels = unique(mu_apatch,'stable');
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = mu_apatch == alpha_levels(aa);
            p = patch(sps(pc),...
                [mu_xpatch(alpha_flags),nan],...
                [mu_ypatch(alpha_flags),nan],0,...
                'edgealpha',alpha_levels(aa),...
                'edgecolor','k',...
                'facecolor','none',...
                'linewidth',1.5);
        end
    else
        p = plot(sps(pc),roi2plot_time,...
            s2_score(:,pc),...
            'color','k',...
            'linestyle','-',...
            'linewidth',1.5);
    end

    % plot projection onset
    onset_flags = roi2plot_time <= 0 & ...
        [roi2plot_time(2:end),nan] > 0;
    plot(sps(pc),roi2plot_time(onset_flags),...
        s2_score(onset_flags,pc),...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',6,...
        'markerfacecolor','w',...
        'markeredgecolor','k');
    
    % iterate through stimuli
    for tt = 1 : n_t
        t2_flags = t2 == t_set(tt);
        trial_flags = ...
            valid_flags & ...
            contrast_flags & ...
            t2_flags;
        if sum(trial_flags) < 1
            continue;
        end
        
        % plot projection offset
        offset_flags = roi2plot_time < t_set(tt) & ...
            [roi2plot_time(2:end),nan] >= t_set(tt);
        scatter(sps(pc),roi2plot_time(offset_flags),...
            s2_score(offset_flags,pc),36,...
            'linewidth',1.5,...
            'marker','o',...
            'markerfacealpha',alpha_levels(tt).^fadeifnoisy,...
            'markerfacecolor','k',...
            'markeredgecolor','none');
    end
end

% update axes
set(sps(pc),...
    'ytick',unique([0,ylim(sps(pc))]),...
    'yticklabel',{'','0',''},...
    'ylim',ylim(sps(pc))+[-1,1]*.15*range(ylim(sps(pc))));

% annotate explained variance
text(sps(pc),...
    .05,.95,sprintf('%.1f%% variance',exp_pca(pc)),...
    'units','normalized');

% ui stacking
uistack(p,'bottom');

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