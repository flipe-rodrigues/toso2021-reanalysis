%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct Ti-aligned, Ii-split psths
pre_padd = 500;
roi2use = [0,t_set(end)];
roi2plot = roi2use; % [-pre_padd,t_set(end)];
roi2use_n_bins = range(roi2use) * psthbin;
roi2plot_n_bins = range(roi2plot) * psthbin;
roi2use_time = linspace(roi2use(1),roi2use(2),roi2use_n_bins);
roi2plot_time = linspace(roi2plot(1),roi2plot(2),roi2plot_n_bins);
roi2use_flags = ...
    roi2plot_time >= roi2use(1) & ...
    roi2plot_time <= roi2use(2);

% preallocation
s1_psths = nan(roi2plot_n_bins,n_neurons);
s2_psths = nan(roi2plot_n_bins,n_neurons);

% distractor flags
i1_flags = i1 == i_set(i1_mode_idx);
i2_flags = i2 == i_set(i2_mode_idx);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    s1_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i1_flags;
    s2_spike_flags = ...
        valid_flags & ...
        neuron_flags & ...
        i2_flags;
    if sum(s1_spike_flags) == 0 || sum(s2_spike_flags) == 0
        continue;
    end
    
    % fetch T1-aligned spike counts & compute spike rates
    s1_spike_counts = data.FR(s1_spike_flags,:);
    s1_spike_rates = conv2(...
        1,kernel.pdf,s1_spike_counts,'valid')' / psthbin * 1e3;
    s1_n_trials = size(s1_spike_counts,1);
    
    % fetch T2-aligned spike counts & compute spike rates
    s2_spike_counts = data.FR(s2_spike_flags,:);
    s2_spike_rates = conv2(...
        1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
    s2_n_trials = size(s2_spike_counts,1);
    
    % T1-aligned spike rates
    s1_alignment_offset = ...
        pre_init_padding + ...
        pre_t1_delay(s1_spike_flags);
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
    
    % T2-aligned spike rates
    s2_alignment = ...
        pre_init_padding + ...
        pre_t1_delay(s2_spike_flags) + ...
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
    s2_psths(:,nn) = nanmean(s2_spkrates,1);
end

%% normalization

% z-score T1-aligned spike density functions
s1_mus = nanmean(s1_psths(roi2use_flags,:),1);
s1_sigs = nanstd(s1_psths(roi2use_flags,:),0,1);
s1_zpsths = (s1_psths - s1_mus) ./ s1_sigs;

% z-score T2-aligned spike density functions
s2_mus = nanmean(s2_psths(roi2use_flags,:),1);
s2_sigs = nanstd(s2_psths(roi2use_flags,:),0,1);
s2_zpsths = (s2_psths - s2_mus) ./ s2_sigs;

%% PCA

% compute T1-aligned PCs
s1_coeff = pca(s1_zpsths(roi2use_flags,:));

% compute T2-aligned PCs
s2_coeff = pca(s2_zpsths(roi2use_flags,:));

%% colormap settings
clim = [-1.5,3];

%% T1-aligned tiling

% figure initialization
fig = figure(figopt,...
    'name','tiling_s1');
axes(axesopt.default,...
    'xlim',roi2plot,...
    'ylim',[1,n_neurons],...
    'xtick',unique([0,t_set',roi2use,roi2plot]),...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8));
title('T1-aligned PSTH raster');
xlabel('Time since T_1 offset (ms)');
ylabel('Neuron #');

% sort by angular position in PC space
[theta,~] = cart2pol(s1_coeff(:,1),s1_coeff(:,2));
[~,theta_idcs] = sortrows(theta);
% theta_idcs = circshift(theta_idcs,sum(theta>0));
theta_idcs = circshift(theta_idcs,80);

% plot selectivity heat map
imagesc(roi2plot,[1,n_neurons],s1_zpsths(:,theta_idcs)',clim);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = ...
    sprintf('z-score_{I_1 = %i %s}',...
    i_set(i1_mode_idx),i1_units);
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

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% T2-aligned tiling

% figure initialization
fig = figure(figopt,...
    'name','tiling_s2');
axes(axesopt.default,...
    'xlim',roi2plot,...
    'ylim',[1,n_neurons],...
    'xtick',unique([0,t_set',roi2use,roi2plot]),...
    'ytick',[1,n_neurons],...
    'colormap',hot(2^8));
title('T2-aligned PSTH raster');
xlabel('Time since T_2 offset (ms)');
ylabel('Neuron #');

% sort by angular position in PC space
[theta,~] = cart2pol(s2_coeff(:,1),s2_coeff(:,2));
[~,theta_idcs] = sortrows(theta);
% theta_idcs = circshift(theta_idcs,sum(theta>0));
theta_idcs = flipud(circshift(theta_idcs,-80));

% plot selectivity heat map
imagesc(roi2plot,[1,n_neurons],s2_zpsths(:,theta_idcs)',clim);

% color bar
clrbar = colorbar;
clrbar.Ticks = unique([0,clim]);
clrlabel.string = ...
    sprintf('z-score_{I_2 = %i %s}',...
    i_set(i2_mode_idx),i2_units);
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

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end