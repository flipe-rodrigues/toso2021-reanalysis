%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% ROI settings
pre_padd = 0;
roi2use = [-pre_padd,isi];
roi2plot = [-pre_padd,isi];
roi2plot_padded = roi2plot + [-1,1] * .05 * range(roi2plot);
roi2use_n_bins = range(roi2use) / psthbin;
roi2plot_n_bins = range(roi2plot_padded) / psthbin;
roi2use_time = linspace(roi2use(1),roi2use(2),roi2use_n_bins);
roi2plot_time = linspace(roi2plot_padded(1),roi2plot_padded(2),roi2plot_n_bins);
roi2use_flags = ...
    roi2plot_time >= roi2use(1) & ...
    roi2plot_time <= roi2use(2);

%% subject selection
subject_flags = ismember(subjects,subject_set);

%% construct isi-aligned psths

% preallocation
zscore_weights = nan(roi2plot_n_bins,n_neurons);
ref_psths = nan(roi2plot_n_bins,n_neurons);
isi_psths = nan(roi2plot_n_bins,n_neurons,n_contrasts);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'parsing neural data');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    ref_spike_flags = ...
        valid_flags & ...
        subject_flags & ...
        neuron_flags;
    if sum(ref_spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    ref_spike_rates = data.SDF(ref_spike_flags,:);
    ref_n_trials = size(ref_spike_rates,1);
    
    % ISI-aligned spike rates
    ref_alignment = ...
        pre_init_padding + ...
        pre_s1_delay(ref_spike_flags) + ...
        t1(ref_spike_flags);
    ref_alignment_flags = ...
        padded_time >= ref_alignment + 0 & ...
        padded_time < ref_alignment + isi;
    ref_chunk_flags = ...
        padded_time >= ref_alignment + roi2plot_padded(1) & ...
        padded_time < ref_alignment + roi2plot_padded(2);
    ref_spkrates = ref_spike_rates';
    ref_spkrates(~ref_alignment_flags') = nan;
    ref_spkrates = reshape(...
        ref_spkrates(ref_chunk_flags'),[roi2plot_n_bins,ref_n_trials])';
    
    % compute observations weights
    zscore_weights(:,nn) = sum(~isnan(ref_spkrates));
    
    % compute mean spike density function
    ref_psths(:,nn) = nanmean(ref_spkrates,1);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        isi_spike_flags = ...
            valid_flags & ...
            subject_flags & ...
            neuron_flags & ...
            contrast_flags;
        
        % fetch spike counts & compute spike rates
        isi_spike_rates = data.SDF(isi_spike_flags,:);
        isi_n_trials = size(isi_spike_rates,1);
        
        % ISI-aligned spike rates
        isi_alignment = ...
            pre_init_padding + ...
            pre_s1_delay(isi_spike_flags) + ...
            t1(isi_spike_flags);
        isi_alignment_flags = ...
            padded_time >= isi_alignment + 0 & ...
            padded_time < isi_alignment + isi;
        isi_chunk_flags = ...
            padded_time >= isi_alignment + roi2plot_padded(1) & ...
            padded_time < isi_alignment + roi2plot_padded(2);
        isi_spkrates = isi_spike_rates';
%         isi_spkrates(~isi_alignment_flags') = nan;
        isi_spkrates = reshape(...
            isi_spkrates(isi_chunk_flags'),[roi2plot_n_bins,isi_n_trials])';
        
        % compute mean spike density function
        isi_psths(:,nn,ii) = nanmean(isi_spkrates,1);
    end
end

% deal with missing entries due to subject selection
nonsubject_flags = all(isnan(ref_psths),1);
ref_psths = ref_psths(:,~nonsubject_flags);
isi_psths = isi_psths(:,~nonsubject_flags,:);
n_selected_neurons = sum(~nonsubject_flags);

% nan handling
if ~strcmpi(contrast_str,'t1')
    isi_psths(isnan(isi_psths)) = 0;
end

%% normalization
mus = nanmean(ref_psths(roi2use_flags,:),1);
sigs = nanstd(ref_psths(roi2use_flags,:),0,1);

% z-scoring
ref_zpsths = (ref_psths - mus) ./ sigs;
isi_zpsths = (isi_psths - mus) ./ sigs;

%% PCA

% training settings
pca_design = ref_zpsths(roi2use_flags,:,:);

% PCA
[coeff,~,~,~,exp_pca] = pca(pca_design);
% coeff = coeff_s2;
ref_score = ref_zpsths * coeff;

% preallocation
isi_score = nan(roi2plot_n_bins,n_selected_neurons,n_contrasts);

% iterate through contrasts
for ii = 1 : n_contrasts
    
    % project onto PCs
    isi_score(:,:,ii) = isi_zpsths(:,:,ii) * coeff;
end

%% recompute explained variance

% preallocation
isi_concat = nan(roi2use_n_bins*n_contrasts,n_selected_neurons);

% concatenate psths across conditions
for nn = 1 : n_selected_neurons
    nn_zpsths_all = isi_zpsths(roi2use_flags,nn,:);
    isi_concat(:,nn) = nn_zpsths_all(:);
end

% compute variance explained
lat_pca = nanvar(isi_concat * coeff)';
exp_pca = lat_pca ./ sum(lat_pca) * 100;

%% 3D trajectories in PC space
fig = figure(figopt,...
    'position',[100,150,460,480],...
    'name',sprintf('pc_trajectories3D_isi_%s',contrast_str));
set(gca,...
    axesopt.default,...
    'xtick',0,...
    'ytick',0,...
    'ztick',0);
xlabel(sprintf('%s\n%.1f%% variance','PC 1',exp_pca(1)),...
    'horizontalalignment','center');
ylabel(sprintf('%s\n%.1f%% variance','PC 2',exp_pca(2)),...
    'horizontalalignment','center');
zlabel(sprintf('%s\n%.1f%% variance','PC 3',exp_pca(3)),...
    'horizontalalignment','center');

% linewidth settings
linewidth = 1.5;
wall_linewidth = linewidth * 2 / 3;

% marker size settings
markersize = 5;
wall_markersize = markersize * 2 / 3;

% iterate through contrasts
for ii = 1 : n_contrasts
    contrast_flags = contrasts == contrast_set(ii);
    
    % plot trajectory
    plot3(isi_score(:,1,ii),...
        isi_score(:,2,ii),...
        isi_score(:,3,ii),...
        'color',contrast_clrs(ii,:),...
        'linestyle','-',...
        'linewidth',linewidth);
    
    % plot S1 offset
    isi_offset_flags = roi2plot_time <= 0 & ...
        [roi2plot_time(2:end),nan] > 0;
    plot3(isi_score(isi_offset_flags,1,ii),...
        isi_score(isi_offset_flags,2,ii),...
        isi_score(isi_offset_flags,3,ii),...
        'linewidth',linewidth,...
        'marker','o',...
        'markersize',markersize,...
        'markerfacecolor','w',...
        'markeredgecolor',contrast_clrs(ii,:));
    
    % plot S2 onset
    s2_onset_flags = roi2plot_time <= isi & ...
        [roi2plot_time(2:end),nan] > isi;
    plot3(isi_score(s2_onset_flags,1,ii),...
        isi_score(s2_onset_flags,2,ii),...
        isi_score(s2_onset_flags,3,ii),...
        'linewidth',linewidth,...
        'marker','o',...
        'markersize',markersize,...
        'markerfacecolor',contrast_clrs(ii,:),...
        'markeredgecolor',contrast_clrs(ii,:));
end

% update axis
axis tight;
spacer = .25;
xlim(xlim + [-1,1] * spacer * range(xlim));
ylim(ylim + [-1,1] * spacer * range(ylim));
zlim(zlim + [-1,1] * spacer * range(zlim));
set(gca,...
    'xtick',[],...xlim,...
    'ytick',[],...ylim,...
    'ztick',[],...zlim,...
    'xticklabel',{},...
    'yticklabel',{},...
    'zticklabel',{},...
    'xcolor','k',...
    'ycolor','k',...
    'zcolor','k');

% plot 2D wall projections
lims = [xlim; ylim; zlim];
% lims(1,:) = fliplr(lims(1,:));
lims(2,:) = fliplr(lims(2,:));
% lims(3,:) = fliplr(lims(3,:));

% plot cube edges
cube_clr = [1,1,1] * .85 * 0;
cube_edges = gobjects(3,1);
cube_edges(1) = plot3([1,1].*lims(1,:),[1,1].*lims(2,1),[1,1].*lims(3,1),...
    'linestyle',':',...
    ...'linewidth',1,...
    'color',cube_clr);
cube_edges(2) = plot3([1,1].*lims(1,1),[1,1].*lims(2,:),[1,1].*lims(3,1),...
    'linestyle',':',...
    ...'linewidth',1,...
    'color',cube_clr);
cube_edges(3) = plot3([1,1].*lims(1,1),[1,1].*lims(2,1),[1,1].*lims(3,:),...
    'linestyle',':',...
    ...'linewidth',1,...
    'color',cube_clr);

% plot pair-wise projections
for jj = 1 : 3
    wall_ref_score = ref_score(:,1:3);
    wall_ref_score(:,jj,:) = lims(jj,1);
    wall_ii_score = isi_score(:,1:3,:);
    wall_ii_score(:,jj,:) = lims(jj,1);
    
    % plot reference trajectory
    %     plot3(wall_ref_score(:,1),...
    %         wall_ref_score(:,2),...
    %         wall_ref_score(:,3),...
    %         'color',[1,1,1]*4/5,...
    %         'linestyle','-',...
    %         'linewidth',1.5);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        
        % faded color
        clrs = colorlerp([contrast_clrs(ii,:);[1,1,1]],5);
        faded_clr = clrs(4,:);
        
        % plot trajectory
        plot3(wall_ii_score(:,1,ii),...
            wall_ii_score(:,2,ii),...
            wall_ii_score(:,3,ii),...
            'color',faded_clr,...
            'linestyle','-',...
            'linewidth',wall_linewidth);
        
        % plot stimulus onset
        isi_offset_flags = roi2plot_time <= 0 & ...
            [roi2plot_time(2:end),nan] > 0;
        plot3(wall_ii_score(isi_offset_flags,1,ii),...
            wall_ii_score(isi_offset_flags,2,ii),...
            wall_ii_score(isi_offset_flags,3,ii),...
            'linewidth',wall_linewidth,...
            'marker','o',...
            'markersize',wall_markersize,...
            'markerfacecolor','w',...
            'markeredgecolor',faded_clr);
        
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
            
            % plot stimulus offset
            s2_onset_flags = roi2plot_time < t_set(tt) & ...
                [roi2plot_time(2:end),nan] >= t_set(tt);
            plot3(wall_ii_score(s2_onset_flags,1,ii),...
                wall_ii_score(s2_onset_flags,2,ii),...
                wall_ii_score(s2_onset_flags,3,ii),...
                'linewidth',wall_linewidth,...
                'marker','o',...
                'markersize',wall_markersize,...
                'markerfacecolor',faded_clr,...
                'markeredgecolor',faded_clr);
            
            % plot cross-condition offset-connecting line
            if ii == contrast_mode_idx
                p = plot3(squeeze(wall_ii_score(s2_onset_flags,1,:)),...
                    squeeze(wall_ii_score(s2_onset_flags,2,:)),...
                    squeeze(wall_ii_score(s2_onset_flags,3,:)),...
                    'linewidth',wall_linewidth*2/3,...
                    'linestyle','-',...
                    'color',faded_clr);
                uistack(p,'bottom');
            end
        end
    end
end

% update axis
view(45,45)

% ui restacking
uistack(cube_edges,'bottom');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% PC projections

% figure initialization
fig = figure(figopt,...
    'position',[535,130,966,860],...
    'name',sprintf('pc_projections_isi_%s',contrast_str));
n_pcs2plot = 6;
sps = gobjects(n_pcs2plot,1);
for pc = 1 : n_pcs2plot
    sp_idx = pc * 2 - 1 - (pc > n_pcs2plot / 2) * (n_pcs2plot - 1);
    sps(pc) = subplot(n_pcs2plot/2,2,sp_idx);
    xlabel(sps(pc),'Time since ISI onset (ms)');
    ylabel(sps(pc),sprintf('PC %i',pc));
end
set(sps,...
    axesopt.default,...
    'xlim',roi2plot + [-1,1] * .05 * range(roi2plot),...
    'ylimspec','tight',...
    'plotboxaspectratio',[3,1,1]);

% link axes
linkaxes(sps,'x');

% iterate through pcs
for pc = 1 : n_pcs2plot
    
    % graphical object preallocation
    p = gobjects(n_contrasts,1);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(ii);
        trial_flags = ...
            valid_flags & ...
            contrast_flags;

        % patch projection
        mu_xpatch = roi2plot_time;
        mu_ypatch = isi_score(:,pc,ii)';
        mu_ypatch(end) = nan;
        p(ii) = plot(sps(pc),roi2plot_time,...
            isi_score(:,pc,ii),...
            'color',contrast_clrs(ii,:),...
            'linestyle','-',...
            'linewidth',1.5);
        
        % plot S1 offset
        isi_offset_flags = roi2plot_time <= 0 & ...
            [roi2plot_time(2:end),nan] > 0;
        plot(sps(pc),roi2plot_time(isi_offset_flags),...
            isi_score(isi_offset_flags,pc,ii),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',6,...
            'markerfacecolor','w',...
            'markeredgecolor',contrast_clrs(ii,:));
        
        % plot S2 onset
        s2_onset_flags = roi2plot_time <= isi & ...
            [roi2plot_time(2:end),nan] > isi;
        plot(sps(pc),roi2plot_time(s2_onset_flags),...
            isi_score(s2_onset_flags,pc,ii),...
            'linewidth',1.5,...
            'marker','o',...
            'markersize',sqrt(36),...
            'markerfacecolor',contrast_clrs(ii,:),...
            'markeredgecolor','none');
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
    uistack(p(:),'bottom');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% 2D trajectories in curated PC space
fig = figure(figopt,...
    'name',sprintf('pc_trajectories2D_isi_%s',contrast_str));
axes(...
    axesopt.default,...
    'clipping','off');
xlabel('Neural dimension 1_{1}');
ylabel('Neural dimension 2_{2}');

% rotation specification
rotfuns = {@rotxd,@rotyd,@rotzd};
thetas = [50,-10,90] .* [1,1,1];
S_ref = nan(3,roi2plot_n_bins,n_contrasts);
S = nan(3,roi2plot_n_bins,n_contrasts);
for pc = 1 : 3
    S_ref(pc,:,1) = ref_score(:,pc);
    S(pc,:,:) = isi_score(:,pc,:);
end
B = cat(3,...
    padarray([-1,0,1],2,0,'post'),...
    padarray([-1,0,1],1,0,'both'),...
    padarray([-1,0,1],2,0,'pre')) * 3;

% preallocate rotation matrix
R = eye(3);

% iterate through dimensions
for ii = 1 : size(R,1)
    
    % compute compound rotation
    R = R * rotfuns{ii}(thetas(ii));
end

% iterate through contrast conditions
for ii = 1 : n_contrasts
    
    % apply rotations
    S_ref(:,:,ii) = R * S_ref(:,:,ii);
    S(:,:,ii) = R * S(:,:,ii);
end

% iterate through bases
for ii = 1 : 3
    
    % apply rotations
    B(:,:,ii) = R * B(:,:,ii);
end

% graphical object preallocation
p = gobjects(n_contrasts,1);

% linewidth settings
linewidth = 1.5;

% marker size settings
markersize = 7.5;

% plot reference trajectory
h_ref = plot3(S_ref(1,:,1),...
    S_ref(2,:,1),...
    S_ref(3,:,1),...
    'color','k',...
    'linestyle','--',...
    'linewidth',linewidth);

% update axis
axis tight;
xxlim = xlim + [-1,1] * .1 * range(xlim);
yylim = ylim + [-1,1] * .1 * range(ylim);
spacer = .05;
xlim(xxlim + [-1,1] * spacer * range(xxlim));
ylim(yylim + [-1,1] * spacer * range(yylim));
set(gca,...
    'xtick',unique([xxlim,0]),...
    'ytick',unique([yylim,0]),...
    'xticklabel',{'','0',''},...
    'yticklabel',{'','0',''},...
    'xcolor','k',...
    'ycolor','k');

% delete reference trajectory
delete(h_ref);

% iterate through contrasts
for ii = 1 : n_contrasts
    contrast_flags = contrasts == contrast_set(ii);
    trial_flags = ...
        valid_flags & ...
        contrast_flags;

    % patch trajectory
    mu_xpatch = S(1,:,ii);
    mu_ypatch = S(2,:,ii);
    mu_ypatch(end) = nan;
    plot(S(1,:,ii),...
        S(2,:,ii),...
        'color',contrast_clrs(ii,:),...
        'linestyle','-',...
        'linewidth',linewidth);
    
    % plot S1 offset
    isi_offset_flags = roi2plot_time <= 0 & ...
        [roi2plot_time(2:end),nan] > 0;
    p(ii) = plot(S(1,isi_offset_flags,ii),...
        S(2,isi_offset_flags,ii),...
        'linewidth',linewidth,...
        'marker','o',...
        'markersize',markersize,...
        'markerfacecolor','w',...
        'markeredgecolor',contrast_clrs(ii,:));
    
    % plot S2 onset
    s2_onset_flags = roi2plot_time <= isi & ...
        [roi2plot_time(2:end),nan] > isi;
    plot(S(1,s2_onset_flags,ii),...
        S(2,s2_onset_flags,ii),...
        'linewidth',1.5,...
        'marker','o',...
        'markersize',markersize,...
        'markerfacecolor',contrast_clrs(ii,:),...
        'markeredgecolor','none');
end

% ui restacking
uistack(p,'top');

% iterate through dimensions
B_clrs = [1,0,0; 0,1,0; 0,0,1];
for ii = 1 : size(B,1)
    
    % plot basis
    plot(B(1,:,ii),...
        B(2,:,ii),...
        'color',B_clrs(ii,:),...
        'linewidth',1.5,...
        'linestyle','-');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end