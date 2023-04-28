%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% run settings
n_runs = 1;

%% selection criteria
n_trial_cutoff = 0;

%% window settings
window_dur = t_set(t2_mode_idx);
t2_flags = t2 >= window_dur;

%% ROI settings

% initialization
roi_win = struct();
roi_lbl = struct();

% roi definitions
roi_win.pre = [-window_dur,0];
roi_win.post = [0,window_dur];

% roi labels
roi_lbl.pre = 'Pre-S2 onset';
roi_lbl.post = 'Post-S2 onset';

% roi length
roi_length.pre = range(roi_win.pre) / psthbin;
roi_length.post = range(roi_win.post) / psthbin;

% parse roi epochs
roi_epochs = fieldnames(roi_win);
n_epochs = numel(roi_epochs);

%% construct spike rate tensor (time X neurons X concatenations)

% preallocation
scale_factors = nan(n_neurons,n_contrasts,n_runs,n_epochs);

% iterate through epochs
for ee = 1 : n_epochs
    epoch = roi_epochs{ee};
    
    % iterate through runs
    for rr = 1 : n_runs
        
        % preallocation
        R = nan(n_neurons,n_contrasts);
        
        % iterate through units
        for nn = 1 : n_neurons
            progressreport(nn,n_neurons,...
                sprintf('fetching spikes (run %i/%i, %s)',rr,n_runs,epoch));
            neuron_flags = data.NeuronNumb == flagged_neurons(nn);
            
            %
            contrast_flags = contrasts ~= -1; %  == contrast_set(contrast_mode_idx);
            
            % trial selection
            trial_flags = ...
                valid_flags & ...
                neuron_flags & ...
                ...t2_flags & ...
                contrast_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials <= n_trial_cutoff
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_counts = data.FR(trial_flags,:);
            spike_rates = ...
                conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
            
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            spike_rates = downsamplecounts(spike_counts,25);
            spike_rates = spike_rates(:,validtime_flags)';
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            % S2-onset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi;
            alignment_flags = ...
                valid_time >= alignment + roi_win.(epoch)(1) & ...
                valid_time < alignment + t2(trial_flags);
            chunk_flags = ...
                valid_time >= alignment + roi_win.(epoch)(1) & ...
                valid_time < alignment + roi_win.(epoch)(2);
            aligned_spkrates = spike_rates;
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [roi_length.(epoch),n_flagged_trials])';
            
            % store tensor & concatenation data
            R(nn,1) = nanmean(aligned_spkrates,[1,2]);

            % iterate through contrasts
            for ii = 1 : n_contrasts
                contrast_flags = contrasts == contrast_set(ii);
                
                % trial selection
                trial_flags = ...
                    valid_flags & ...
                    neuron_flags & ...
                    ...t2_flags & ...
                    contrast_flags;
                flagged_trials = find(trial_flags);
                n_flagged_trials = numel(flagged_trials);
                if n_flagged_trials <= n_trial_cutoff
                    continue;
                end
                
                % fetch spike counts & compute spike rates
                spike_counts = data.FR(trial_flags,:);
                spike_rates = ...
                    conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
                
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        spike_rates = downsamplecounts(spike_counts,25);
        spike_rates = spike_rates(:,validtime_flags)';
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
                % S2-onset-aligned spike rates
                alignment = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags) + ...
                    isi;
                alignment_flags = ...
                    valid_time >= alignment + roi_win.(epoch)(1) & ...
                    valid_time < alignment + t2(trial_flags);
                chunk_flags = ...
                    valid_time >= alignment + roi_win.(epoch)(1) & ...
                    valid_time < alignment + roi_win.(epoch)(2);
                aligned_spkrates = spike_rates;
                aligned_spkrates(~alignment_flags') = nan;
                aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                    [roi_length.(epoch),n_flagged_trials])';

                % store tensor & concatenation data
                R(nn,ii+1) = nanmean(aligned_spkrates,[1,2]);
            end
        end
        
        % compute scale factors
        scale_factors(:,:,rr,ee) = R(:,2:end) ./ R(:,1);
    end
    
    %% plot response dilation
    
    % set dilation boundaries
    dilation_bounds = [-1,1] * 30;
    
    % figure initialization
    fig = figure(figopt,...
        ...'position',[250,730,420,315],...
        'name',sprintf('avgfr_dilation_%s_%s',epoch,contrast_str));
    axes(axesopt.default,...
        'xlim',[1,n_contrasts] + [-1,1]*.65,...
        'xtick',1:n_contrasts,...
        'xticklabel',num2cell(contrast_set),...
        'ylim',dilation_bounds+[-1,1]*range(dilation_bounds)*.05,...
        'ytick',linspace(dilation_bounds(1),dilation_bounds(2),5),...
        'xticklabelrotation',45,...
        ...'ticklength',axesopt.default.ticklength/.75,...
        'layer','bottom',...
        'plotboxaspectratio',[1,1.5,1]);
    title(capitalize(epoch));
    xlabel(sprintf('%s (%s)',contrast_lbl,contrast_units));
    ylabel('<Firing rate> dilation (%)');
    
    % convert scaling factors to dilations
    all_dilations = (nanmean(scale_factors(:,:,:,ee),3) - 1) * 100;
%     all_dilations = nanmean(scale_factors(:,:,:,ee),3) * 100;
%     all_dilations = R(:,2:end) * 1.5e0;
%     all_dilations = squeeze(nanmean(zpsths(time>=-334&time<=0,:,:),1)) * 30;
%     all_dilations = squeeze(nanmean(zpsths(time>=0&time<=334,:,:),1)) * 30;
    ref_dilations = all_dilations(:,contrast_mode_idx);
    
    % zero line
    plot(xlim,[1,1]*0,...
        'color','k',...
        'linestyle',':');
    
    % graphical object preallocation
    p = gobjects(3,n_contrasts);
    
    % bar settings
    barwidth = 2/3;
    
    % dilation selection settings
    lowerbound = min(ylim);
    upperbound = max(ylim);
    % lowerbound = -inf;
    % upperbound = +inf;
    
    % average function selection
    avgfun = @(x,d) nanmedian(x,d);
    errfun = @(x) diff(quantile(x,3));
%     avgfun = @(x,d) nanmean(x,d);
%     errfun = @(x) [1,1] * nanstd(x) / sqrt(numel(x));
    
    % compute the reference average
    fov_flags = all(...
        all_dilations >= lowerbound & ...
        all_dilations <= upperbound,2);
    offset = avgfun(ref_dilations,1); % avgfun(dilations(fov_flags,:),[1,2]);

    %
    dilation_flags = ...
        ...fov_flags & ...
        true(n_neurons,1);
    n_flaggeddilations = sum(dilation_flags)
    dilation_edges = linspace(...
        dilation_bounds(1)*1.25,dilation_bounds(2)*1.25,n_flaggeddilations);
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        
        % mean & standard error across animals
        cond_avg = avgfun(all_dilations(dilation_flags,ii),1) - offset;
        cond_err = errfun(all_dilations(dilation_flags,ii));
        
        % plot summary of threshold dilation
        xpatch = ii + [-1,1,1,-1] * barwidth / 2;
        ypatch = cond_avg + [[-1,-1] * cond_err(1),[+1,+1] * cond_err(2)];
        
        p(1,ii) = plot(ii+[-1,1]*.65*barwidth,[1,1]*cond_avg,...
            'color',contrast_clrs(ii,:),...
            'linestyle','-',...
            'linewidth',2.5);
        p(2,ii) = plot(ii+[-1,1]*.5*barwidth,[1,1]*cond_avg,...
            'color','w',...
            'linestyle','-',...
            'linewidth',5);
        p(3,ii) = patch(xpatch,ypatch,contrast_clrs(ii,:),...
            'edgecolor','w',...
            'facealpha',1,...
            'linewidth',1);
        
        % plot ground truth
        plot(xlim,[1,1]*(contrast_set(ii)-1)*100,...
            'color',contrast_clrs(ii,:),...
            'linestyle',':',...
            'linewidth',1);
    end
    
    % iterate through contrasts
    for ii = 1 : n_contrasts
        [dilation_pdf,~] = ksdensity(...
            all_dilations(dilation_flags,ii) - offset,dilation_edges);
        dilation_pdf = 4/15 * barwidth * ...
            (dilation_pdf - min(dilation_pdf)) ./ range(dilation_pdf);
        noise = randn(n_flaggeddilations,1) .* dilation_pdf';
        
        % plot animal threshold dilation
        grapeplot(ii+noise,...
            sort(all_dilations(dilation_flags,ii)) - offset,...
            'marker','o',...
            'markersize',4,...
            'markeredgecolor',contrast_clrs(ii,:),...
            'markerfacecolor',[1,1,1],...
            'linewidth',1);
    end
    
    % ui restacking
    sp = p(1:2,:);
    uistack(sp(:),'top');
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end

    %% compute response stretch
    
    % preallocation
    stretch = nan(n_neurons,1);
    
    % design matrix
    x = contrast_set;
    X = [ones(n_contrasts,1),x];
    
    % iterate through units
    for uu = 1 : n_neurons
        
        % linear regression
        thetas = X \ (all_dilations(uu,:)' - offset);
        stretch(uu) = thetas(end);
    end
    
    % one-sample t-test
    [~,pvalue,~,stats.ttest] = ttest(stretch);
    stats.ttest.pvalue = pvalue;
    
    %% plot response stretch
    stretch_bounds = [-1,1] * 1 / (10^strcmpi(contrast_str,'t1'));
    count_bounds = [0,75];
    
    % figure initialization
    fig = figure(figopt,...
        'name',sprintf('avgfr_stretch_%s_%s',epoch,contrast_str));
    axes(...
        axesopt.default,...
        'plotboxaspectratio',[1,3,1],...
        'xlim',count_bounds,...
        'xtick',count_bounds,...
        'xdir','reverse',...
        'ylim',stretch_bounds+[-1,1]*.05*range(stretch_bounds),...
        'ytick',linspace(stretch_bounds(1),stretch_bounds(2),5),...
        'yaxislocation','right');
    title(capitalize(epoch));
    xlabel('Count');
    ylabel(sprintf('Response stretch (%% / %s)',contrast_units),...
        'rotation',-90,...
        'verticalalignment','bottom');
    
    % zero line
    plot(xlim,[1,1]*0,...
        'color','k',...
        'linestyle',':');
    
    % distribution settings
    nbins = 41;
    bin_edges = linspace(min(yticks),max(yticks),nbins);
    bin_counts = histcounts(stretch,bin_edges);
    
    % plot stretch distribution
    histogram(...
        'binedges',bin_edges,...
        'bincounts',bin_counts,...
        'orientation','horizontal',...
        'edgecolor','k',...
        'facecolor','k',...
        'facealpha',.5,...
        'linestyle','none');
    stairs([0,bin_counts],bin_edges,...
        'color','k',...
        'linewidth',1.5);
    
    % compute mean stretch
    fov_flags = ...
        stretch >= min(ylim) & ...
        stretch <= max(ylim);
    stretch_avg = nanmean(stretch(fov_flags));
    
    % assess population significance
    if pvalue <= .01
        str = '**';
    elseif pvalue <= .05
        str = '*';
    else
        str = '^{n.s.}';
    end
    text(max(xlim)*.9,.025*range(ylim)+stretch_avg,str,...
        'color','k',...
        'fontname',axesopt.default.fontname,....
        'fontsize',axesopt.default.fontsize*1.25,...
        'horizontalalignment','center');
    
    % test affordance
    p = plot([1,1]*max(xlim)*.9,[0,1]*stretch_avg,...
        'color','k',...
        'linewidth',1.5);
    uistack(p,'bottom');
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end