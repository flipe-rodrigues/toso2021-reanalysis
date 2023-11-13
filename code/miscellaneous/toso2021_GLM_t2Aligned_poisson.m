%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% GLM settings
distro = 'poisson';
glm_win = 250;
glm_roi = [-500-glm_win,t_set(end)];
glm_step = 25;
n_glm = floor((diff(glm_roi) - glm_win) / glm_step) + 1;
glm_time = linspace(glm_roi(1)+glm_win,glm_roi(2)-glm_step,n_glm);

%% neuron selection

% selected for being good examples of i2-modulation
if strcmpi(task_str,'duration')
    neurons2use = [...
        21,24,35,38,62,65,68,72,100,130,205,206,215,224,...
        234,241,356,381,391,393,397,402,406,428,441,448,...
        459,461,462,470,473,493,526,544,553,555,566];
%     neurons2use = [...
%         38,72,205,215,224,391,393,397,402,448,459,462,470,526,566];
elseif strcmpi(task_str,'intensity')
    neurons2use = [...
        19,22,30,61,66,70,100,111,112,115,...
        166,238,243,260,344,408,410];
end
neurons2use = flagged_neurons;
neurons2use = neuron_idcs;
n_neurons2use = numel(neurons2use);

%% construct response

% preallocation
spkcounts_s2 = nan(n_total_trials,n_glm);

% iterate through neurons
for nn = 1 : n_neurons2use
    progressreport(nn,n_neurons2use,'fetching spike counts');
    neuron_flags = data.NeuronNumb == neurons2use(nn);
    
    % flag trials for the current condition
    spike_flags = ...
        valid_flags & ...
        neuron_flags;
    flagged_trials = find(spike_flags);
    if sum(spike_flags) == 0
        continue;
    end
    
    % fetch spike counts & compute spike rates
    spike_counts = data.FakeFR(spike_flags,:)';
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
            pre_s1_delay(spike_flags) + ...
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
        spkcounts_s2(spike_flags,gg) = nansum(spkrates,2);
    end
end

%% assess glm validity through time

% preallocation
valid_glms = false(n_neurons2use,n_glm);

% iterate through neurons
for nn = 1 : n_neurons2use
    progressreport(nn,n_neurons2use,'asessing validity');
    
    % fetch current neuron index
    nn_idx = neurons2use(nn);
    
    % iterate through GLMs
    for gg = 1 : n_glm
        glm_win_onset = glm_roi(1) + (gg - 1) * glm_step;
        glm_win_offset = glm_win_onset + glm_win;
        
        % fetch current stimulus index
        t2_idx = find(glm_win_offset <= t_set,1);
        
        % assess validity
        valid_glms(nn,gg) = ...
            all(surviving_trial_counts(nn_idx,t2_idx,:) > trial_count_cutoff);
    end
end

%% spike count GLMs

% design matrix
X = [s1,d1,d2,trial_idcs];
% X = [s1,d1,d2];
X = [d1,d2];
% X = [d2];
n_regressors = size(X,2);
n_coefficients = n_regressors + 1;
    
% preallocation
betas_s2 = zeros(n_neurons2use,n_glm,n_coefficients);
pvals_s2 = zeros(n_neurons2use,n_glm,n_coefficients);

% feature normalization
Z = (X - nanmean(X)) ./ nanstd(X);

% iterate through neurons
for nn = 1 : n_neurons2use
    progressreport(nn,n_neurons2use,'fitting neuron-wise GLMs');
    neuron_flags = data.NeuronNumb == neurons2use(nn);
    trial_flags = ...
        valid_flags & ...
        neuron_flags;
    if sum(trial_flags) == 0
        continue;
    end

    % iterate through GLMs
    for gg = 1 : n_glm

        % fit GLM to each subject
        mdl = fitglm(Z(trial_flags,:),spkcounts_s2(trial_flags,gg),'linear',...
            ...'predictorvars',{s1_lbl,d1_lbl,d2_lbl,'trial #'},...
            ...'predictorvars',{s1_lbl,d1_lbl,d2_lbl,'trial #'},...
            'predictorvars',{d1_lbl,d2_lbl},...
            'distribution',distro,...
            'intercept',true);
        betas_s2(nn,gg,:) = mdl.Coefficients.Estimate;
        pvals_s2(nn,gg,:) = mdl.Coefficients.pValue;
    end
end

%% print percentage of significantly modulated neurons
alpha = .01;
significance_mask = squeeze(pvals_s2(:,:,2:end)) < alpha;
n_significant = sum(squeeze(nansum(significance_mask,2)) >= 1);

% iterate through coefficients
for bb = 1 : n_regressors
    fprintf('%s-modulated: %i/%i (%.1f%%)\n',...
        mdl.Coefficients.Properties.RowNames{bb+1},...
        n_significant(bb),n_neurons2use,...
        n_significant(bb)/n_neurons2use*100);
end

%% color settings
r_clr = [.95, .25, .25];
w_clr = [1, 1, 1] * 1;
b_clr = [.15, .35, .75];
clrmap = colorlerp([b_clr; w_clr; r_clr], 2^8);

%% neuron highlights
eg_neurons = [68,72,215,391,393,428,459,470,526];
n_egneurons = numel(eg_neurons);

%% GLM significance heatmaps

% colormap settings
clims = [-2,3];

% iterate through coefficients
for bb = 2 : n_coefficients
    coeff_lbl = mdl.Coefficients.Properties.RowNames{bb};
    coeff_flags = ismember(mdl.Coefficients.Properties.RowNames,coeff_lbl);
    coeff_lbl = strrep(coeff_lbl,' #','');

    % coefficient map
    coeff_map_s2.(coeff_lbl) = betas_s2(:,:,coeff_flags)';
    
    % significance map
    significance_mask_s2.(coeff_lbl) = ...
        double(pvals_s2(:,:,coeff_flags)' < .05) .* sign(coeff_map_s2.(coeff_lbl)) + ...
        double(pvals_s2(:,:,coeff_flags)' < .01) .* sign(coeff_map_s2.(coeff_lbl));
%     significance_mask(~valid_glms') = max(clims);
    
    % pca
    [pc_coeff,~,~,~,exp] = pca(significance_mask_s2.(coeff_lbl));
    exp_cutoff = 80;
    n_pcs2use = max(sum(cumsum(exp) <= exp_cutoff),2);
    
    % average sorting
    avg_coeffs = mean(coeff_map_s2.(coeff_lbl));
    [~,avg_idcs] = sort(avg_coeffs);
    
    % sort by angular position in PC space
    [theta,~] = cart2pol(pc_coeff(:,1),pc_coeff(:,2));
    [~,theta_idcs] = sortrows(theta);
    
    % agglomerative hierarchical clustering
%     diss = pdist(pc_coeff(:,1:n_pcs2use),'euclidean');
    diss = pdist(significance_mask_s2.(coeff_lbl)','euclidean');
    tree = linkage(diss,'ward');
    leaf_idcs = optimalleaforder(tree,diss);

    % color map settings
    coeff_str = strrep(lower(coeff_lbl),'_','');
    coeff_clrs = eval([coeff_str,'_clrs']);
    clrmap = [...
        colorlerp([coeff_clrs(1,:); [1,1,1]; coeff_clrs(end,:)], 5);...
        [0,0,1]];
    
    % figure initialization
    fig = figure(figopt,...
        'position',[1.8+(bb-2)*405,41.8,405,1.0288e+03],...
        'name',sprintf('GLM_significance_%s',...
        strrep(lower(coeff_lbl),'_','')));
    xxtick = unique([glm_roi';glm_roi(1)+glm_win;0;t_set]);
    xxticklabel = num2cell(xxtick);
    xxticklabel(~ismember(xxtick,[glm_roi,glm_roi(1)+glm_win,0])) = {''};
    axes(axesopt.default,...
        'plotboxaspectratiomode','auto',...
        'ticklength',[0.0250 0.0250]/3,...
        'xlim',glm_roi+[1,0]*glm_win,...
        'ylim',[1,n_neurons2use],...
        'xtick',xxtick,...
        'xticklabel',xxticklabel,...
        'ytick',[1,n_neurons2use],...
        'layer','top',...
        'clipping','off',...
        'colormap',clrmap);
%     title(coeff_lbl);
    xlabel('Time since T_2 onset (ms)');
    ylabel('Neuron #');

    % plot selectivity heat map
    sorted_idcs = leaf_idcs;
%     significance_mask(valid_glms') = 0;
%     significance_mask(~valid_glms') = max(clims);
    significance_mask_s2.(coeff_lbl) = ...
        medfilt2(significance_mask_s2.(coeff_lbl),[3,1]);
    imagesc(glm_roi+[1,0]*glm_win+[1,-1]*glm_step/2,[1,n_neurons2use],...
        significance_mask_s2.(coeff_lbl)(:,sorted_idcs)',clims);
 
    % plot alignment line
    plot([1,1]*0,ylim,...
        'color','k',...
        'linestyle','--',...
        'linewidth',1);
        
    % iterate through example neurons
%     for ee = 1 : n_egneurons
%         try
%         eg_neuron_flags = neurons2use(sorted_idcs) == eg_neurons(ee);
%         
%         % highlight example neuron
%         plot(min(xlim)-range(xlim)*.0275,find(eg_neuron_flags),...
%             'marker','>',...
%             'markersize',10,...
%             'markeredgecolor','k',...
%             'markerfacecolor',[1,1,1]*(ee-1)/(n_egneurons-1),...
%             'linewidth',1.5);
%         catch
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
    delete(clrbar);
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end