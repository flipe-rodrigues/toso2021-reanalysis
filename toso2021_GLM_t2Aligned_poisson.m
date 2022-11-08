%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% GLM settings
distro = 'poisson';
glm_roi = [-500,t_set(t2_mode_idx+2)];
glm_win = 250;
glm_step = 25;
n_glm = floor((diff(glm_roi) - glm_win) / glm_step) + 1;
glm_time = linspace(glm_roi(1),glm_roi(2)-glm_step,n_glm);

%% construct response

% preallocation
spkcounts = nan(n_total_trials,n_glm);

% clamping
i1_flags = i1 == i_set(i1_mode_idx);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'fetching spike counts');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    
    % flag trials for the current condition
    spike_flags = ...
        valid_flags & ...
        ...i1_flags & ...
        neuron_flags;
    flagged_trials = find(spike_flags);
    if sum(spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    spike_counts = data.FR(spike_flags,:)';
    spike_rates = ...
        conv2(kernel.pdf,1,spike_counts,'valid') / psthbin * 1e3;
    n_trials = sum(spike_flags);
    
    % iterate through GLMs
    for gg = 1 : n_glm
        glm_win_onset = glm_roi(1) + (gg - 1) * glm_step;
        glm_win_offset = glm_win_onset + glm_win;
        
        % T2-onset-aligned spike rates
        alignment_onset = ...
            pre_init_padding + ...
            pre_t1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi;
        alignment_flags = ...
            padded_time >= alignment_onset + glm_roi(1) & ...
            padded_time < alignment_onset + t2(spike_flags);
        chunk_flags = ...
            padded_time >= alignment_onset + glm_win_onset & ...
            padded_time < alignment_onset + glm_win_offset;
        spkrates = spike_counts;
        spkrates(~alignment_flags') = nan;
        spkrates = reshape(spkrates(chunk_flags'),[glm_win,n_trials])';

        % store average spike rates
        spkcounts(spike_flags,gg) = nansum(spkrates,2);
    end
end

%% spike count GLMs

% design matrix
X = [s1,d1,d2,data.Trial];
X = [s1,d1,d2];
% X = [d1,d2];
% X = [d2];
n_regressors = size(X,2);
n_coefficients = n_regressors + 1;
    
% preallocation
betas = nan(n_neurons,n_glm,n_coefficients);
pvals = nan(n_neurons,n_glm,n_coefficients);

% feature normalization
Z = (X - nanmean(X)) ./ nanstd(X);

% iterate through neurons
for nn = 1 : n_neurons
    progressreport(nn,n_neurons,'fitting neuron-wise GLMs');
    neuron_flags = data.NeuronNumb == flagged_neurons(nn);
    trial_flags = ...
        valid_flags & ...
        neuron_flags;
    
    % iterate through GLMs
    for gg = 1 : n_glm

        % fit GLM to each subject
        mdl = fitglm(Z(trial_flags,:),spkcounts(trial_flags,gg),'linear',...
            'predictorvars',{s1_lbl,d1_lbl,d2_lbl},...{s1_lbl,d1_lbl,d2_lbl,'trial #'},...
            'distribution',distro,...
            'intercept',true);
        betas(nn,gg,:) = mdl.Coefficients.Estimate;
        pvals(nn,gg,:) = mdl.Coefficients.pValue;
    end
end

%% print percentage of significantly modulated neurons
alpha = .01;
significance_mask = squeeze(pvals(:,:,2:end)) < alpha;
n_significant = sum(squeeze(nansum(significance_mask(:,glm_time>=0,:),2)) >= 1);

% iterate through coefficients
for bb = 1 : n_regressors
    fprintf('%s-modulated: %i/%i (%.1f%%)\n',...
        mdl.Coefficients.Properties.RowNames{bb+1},...
        n_significant(bb),n_neurons,...
        n_significant(bb)/n_neurons*100);
end

%% color settings
r_clr = [.95, .25, .25];
w_clr = [1, 1, 1] * 1;
b_clr = [.15, .35, .75];
clrmap = colorlerp([b_clr; w_clr; r_clr], 2^8);

%% neuron highlights
eg_neurons = [215,393,526];
eg_neurons = [68,72,215,391,393,428,459,470,526];
n_egneurons = numel(eg_neurons);

%% GLM coefficient heatmaps

% colormap settings
clims = [-1,1] * .75;

% iterate through coefficients
for bb = 2 : n_coefficients
    coeff_lbl = mdl.Coefficients.Properties.RowNames{bb};
    coeff_flags = ismember(mdl.Coefficients.Properties.RowNames,coeff_lbl);

    % coefficient map
    coeff_map = betas(:,:,coeff_flags)';
    
    % coefficient significance
    significance_mask = pvals(:,:,coeff_flags)' < alpha;
  
    % pca
    [pc_coeff,~,~,~,exp] = pca(coeff_map);
    exp_cutoff = 80;
    n_pcs2use = max(sum(cumsum(exp) <= exp_cutoff),2);
    
    % average sorting
    avg_coeffs = mean(coeff_map);
    [~,avg_idcs] = sort(avg_coeffs);
    
    % sort by angular position in PC space
    [theta,~] = cart2pol(pc_coeff(:,1),pc_coeff(:,2));
    [~,theta_idcs] = sortrows(theta);
    
    % agglomerative hierarchical clustering
    diss = pdist(pc_coeff(:,1:n_pcs2use),'euclidean');
    diss = pdist(coeff_map','euclidean');
    tree = linkage(diss,'ward');
    leaf_idcs = optimalleaforder(tree,diss);

    % color map settings
    coeff_str = strrep(lower(coeff_lbl),'_','');
    coeff_clrs = eval([coeff_str,'_clrs']);
    clrmap = colorlerp(...
        [coeff_clrs(1,:); [1,1,1]; coeff_clrs(end,:)], 2^8);
    
    % figure initialization
    fig = figure(figopt,...
        ...'windowstyle','docked',...
        ...'position',[545+(bb-2)*500,912,500,1310],...
        'name',sprintf('GLM_coefficients_%s',coeff_str));
    axes(axesopt.default,...
        'plotboxaspectratiomode','auto',...
        'xlim',glm_roi,...
        'ylim',[1,n_neurons],...
        'xtick',unique([0,t_set',glm_roi]),...
        'ytick',[1,n_neurons],...
        'layer','top',...
        'clipping','off',...
        'colormap',clrmap);
    title(coeff_lbl);
    xlabel('Time since T_2 onset (ms)');
    ylabel('Neuron #');

    % plot selectivity heat map
    sorted_idcs = leaf_idcs;
    coeff_map(~significance_mask) = 0;
%     coeff_map = medfilt2(coeff_map,[3,1]);
    imagesc(glm_roi+[1,-1]*glm_step,[1,n_neurons],...
        coeff_map(:,sorted_idcs)',clims);

    % plot alignment line
    plot([1,1]*0,ylim,...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);

    % iterate through example neurons
    for ee = 1 : n_egneurons
        eg_neuron_flags = flagged_neurons(sorted_idcs) == eg_neurons(ee);
        
        % highlight example neuron
        plot(min(xlim)-range(xlim)*.0275,find(eg_neuron_flags),...
            'marker','>',...
            'markersize',10,...
            'markeredgecolor','k',...
            'markerfacecolor',[1,1,1]*(ee-1)/(n_egneurons-1),...
            'linewidth',1.5);
    end
    
    % color bar
    clrbar = colorbar;
    clrbar.Ticks = unique([0,clims]);
    clrlabel.string = 'Regression weight';
    clrlabel.fontsize = axesopt.default.fontsize * 1.1;
    clrlabel.rotation = 270;
    clrlabel.position = [4.4,0,0];
    clrlabel.position = [4.4,sum(clims)/2,0];
    set(clrbar,...
        axesopt.colorbar,...
        'fontsize',axesopt.default.fontsize);
    set(clrbar.Label,...
        clrlabel);

    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end
return;

%% GLM significance heatmaps

% colormap settings
clims = [-1,1];

% iterate through coefficients
for bb = 2 : n_coefficients
    coeff_lbl = mdl.Coefficients.Properties.RowNames{bb};
    coeff_flags = ismember(mdl.Coefficients.Properties.RowNames,coeff_lbl);

    % coefficient map
    coeff_map = betas(:,:,coeff_flags)';
    
    % significance map
    significance_mask = pvals(:,:,coeff_flags)' < alpha;
    
    % pca
    [pc_coeff,~,~,~,exp] = pca(coeff_map);
    exp_cutoff = 80;
    n_pcs2use = max(sum(cumsum(exp) <= exp_cutoff),2);
    
    % average sorting
    avg_coeffs = mean(coeff_map);
    [~,avg_idcs] = sort(avg_coeffs);
    
    % sort by angular position in PC space
    [theta,~] = cart2pol(pc_coeff(:,1),pc_coeff(:,2));
    [~,theta_idcs] = sortrows(theta);
    
    % agglomerative hierarchical clustering
    diss = pdist(pc_coeff(:,1:n_pcs2use),'euclidean');
    diss = pdist((double(significance_mask) .* sign(coeff_map))','euclidean');
    tree = linkage(diss,'ward');
    leaf_idcs = optimalleaforder(tree,diss);

    % color map settings
    coeff_str = strrep(lower(coeff_lbl),'_','');
    coeff_clrs = eval([coeff_str,'_clrs']);
    clrmap = colorlerp(...
        [coeff_clrs(1,:); [1,1,1]; coeff_clrs(end,:)], 3);
    
    % figure initialization
    fig = figure(figopt,...
        ...'windowstyle','docked',...
        ...'position',[545+(bb-2)*500,912,500,1310],...
        'name',sprintf('GLM_significance_%s',...
        strrep(lower(coeff_lbl),'_','')));
    axes(axesopt.default,...
        'plotboxaspectratiomode','auto',...
        'ticklength',[0.0250 0.0250]/3,...
        'xlim',glm_roi,...
        'ylim',[1,n_neurons],...
        'xtick',unique([0,t_set',glm_roi]),...
        'ytick',[1,n_neurons],...
        'layer','top',...
        'clipping','off',...
        'colormap',clrmap);
    title(coeff_lbl);
    xlabel('Time since S_2 onset (ms)');
    ylabel('Neuron #');

    % plot selectivity heat map
    sorted_idcs = leaf_idcs;
    significance_mask = double(significance_mask) .* sign(coeff_map);
%     significance_mask = medfilt2(significance_mask,[3,1]);
    imagesc(glm_roi,[1,n_neurons],...
        significance_mask(:,sorted_idcs)',clims);
 
    % plot alignment line
    plot([1,1]*0,ylim,...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
        
    % iterate through example neurons
%     for ee = 1 : n_egneurons
%         eg_neuron_flags = flagged_neurons(sorted_idcs) == eg_neurons(ee);
%         
%         % highlight example neuron
%         plot(min(xlim)-range(xlim)*.0275,find(eg_neuron_flags),...
%             'marker','>',...
%             'markersize',10,...
%             'markeredgecolor','k',...
%             'markerfacecolor',[1,1,1]*(ee-1)/(n_egneurons-1),...
%             'linewidth',1.5);
%         if bb == 4
%             a=1
%         end
%     end
    
    % color bar
    clrbar = colorbar;
    clrbar.Ticks = unique([0,clims]);
    clrlabel.string = 'Regression weight';
    clrlabel.fontsize = axesopt.default.fontsize * 1.1;
    clrlabel.rotation = 270;
    clrlabel.position = [4.4,0,0];
    clrlabel.position = [4.4,sum(clims)/2,0];
    set(clrbar,...
        axesopt.colorbar,...
        'fontsize',axesopt.default.fontsize);
    set(clrbar.Label,...
        clrlabel);
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end