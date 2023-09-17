%% initialization
if ~exist('data','var')
    toso2021_wrapper;
    close all;
end

%% bootstrap settings
n_boots = 10;

%% temporal smoothing kernel
gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);

%% time settings
dt = 75;
ti = 0 + gauss_kernel.paddx(1);
tf = t_set(t2_mode_idx+1) + gauss_kernel.paddx(2);
T = (tf - ti) / psthbin;
t = linspace(ti,tf,T);
% t = ti : tf;
% T = numel(t);
t_padded = ti - dt/2 : tf + dt/2;
t_units = 1e3;

%% epoch settings
epochs = {'s1','s2'};
n_epochs = numel(epochs);

%% cluster settings
cluster_labels = {'ramp','nonramp','corr','noncorr'};
n_clusters = numel(cluster_labels);
min_cluster_size = inf;

% preallocation
cluster_idcs = cell(n_epochs,n_clusters);

% iterate through epochs
for ee = 1 : n_epochs
    epoch = epochs{ee};

    % iterate through clusters
    for kk = 1 : n_clusters
        cluster = cluster_labels{kk};
        all_idcs = eval([cluster,'_idcs']).(epoch);
        cluster_idcs{ee,kk} = all_idcs(ismember(all_idcs,flagged_neurons));
        min_cluster_size = min(min_cluster_size,numel(cluster_idcs{ee,kk}));
    end
end

% table conversion
cluster_idcs = cell2table(cluster_idcs',...
    'variablenames',epochs,...
    'rownames',cluster_labels);

% color settings
cluster_clrs = [...
    ramp_clrs;...
    corr_clrs];

%% cross-validation settings
n_folds = 2;

%% construct spike rate tensor (time X neurons X conditions)

% data clearance
clear R P_tR P_tR;

% preallocation
P_tR = nan(T,T,n_epochs,n_clusters,n_boots);
P_tR_mu = nan(T,n_epochs,n_clusters,n_boots);
P_tR_sd = nan(T,n_epochs,n_clusters,n_boots);

% iterate through bootstrap iterations
for bb = 1 : n_boots
    
    % preallocation
    R_all = nan(T,n_neurons_total,n_epochs,n_folds);
    
    % iterate through units
    for nn = 1 : n_neurons_total
        progressreport(nn,n_neurons_total,...
            sprintf('sampling concatenations (boot %i/%i)',bb,n_boots));
        neuron_flags = data.NeuronNumb == neuron_idcs(nn);
        spike_flags = ...
            valid_flags & ...
            neuron_flags;
        if sum(spike_flags) == 0 || ~ismember(nn,flagged_neurons)
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_rates = data.SDF(spike_flags,:);
        n_trials = size(spike_rates,1);
        
        % S1-onset-aligned spike rates
        s1_alignment = ...
            pre_init_padding + ...
            pre_s1_delay(spike_flags);
        s1_alignment_flags = ...
            padded_time >= s1_alignment + ti & ...
            padded_time < s1_alignment + t1(spike_flags);
        s1_chunk_flags = ...
            padded_time >= s1_alignment + ti & ...
            padded_time < s1_alignment + tf;
        s1_spkrates = spike_rates';
        s1_spkrates(~s1_alignment_flags') = nan;
        s1_spkrates = reshape(s1_spkrates(s1_chunk_flags'),...
            [T,n_trials])';
        
        % S2-onset-aligned spike rates
        s2_alignment = ...
            pre_init_padding + ...
            pre_s1_delay(spike_flags) + ...
            t1(spike_flags) + ...
            isi;
        s2_alignment_flags = ...
            padded_time >= s2_alignment + ti & ...
            padded_time < s2_alignment + t2(spike_flags);
        s2_chunk_flags = ...
            padded_time >= s2_alignment + ti & ...
            padded_time < s2_alignment + tf;
        s2_spkrates = spike_rates';
        s2_spkrates(~s2_alignment_flags') = nan;
        s2_spkrates = reshape(s2_spkrates(s2_chunk_flags'),...
            [T,n_trials])';
        
        % handle cross-validation
        train_flags = ismember(1:n_trials,...
            randperm(n_trials,floor(n_trials/n_folds)));
        
        % iterate through epochs
        for ee = 1 : n_epochs
            epoch_spkrates = eval([epochs{ee},'_spkrates']);
            
            % iterate through train & test sets
            for tt = [0,1]
                r = epoch_spkrates(train_flags==tt,:);

                % compute average firing rates through time
                r_mu = nanmean(r);
                nan_flags = isnan(r_mu);

                % temporal smoothing
                r_gauss = conv2(1,gauss_kernel.pdf,r_mu,'valid');
                r_gauss = padarray(r_gauss,[0,floor(gauss_kernel.nbins/2)],nan,'pre');
                r_gauss = padarray(r_gauss,[0,ceil(gauss_kernel.nbins/2)-1],nan,'post');
                r_gauss(nan_flags) = nan;

                % store train & test spike rates
                R_all(:,nn,ee,2-tt) = r_gauss;
            end
        end
    end
    
    %% decode time from each epoch & cluster
    
    % iterate through epochs
    for ee = 1 : n_epochs
        epoch = epochs{ee};
        
        % iterate through clusters
        for kk = 1 : n_clusters
            cluster_size = numel(cluster_idcs.(epoch){kk});
            cluster_flags = ismember(1:cluster_size,...
                randperm(cluster_size,min_cluster_size));
            R = R_all(:,cluster_idcs.(epoch){kk}(cluster_flags),ee,:);
            
            %% naive bayes decoder
            nbdopt = struct();
            nbdopt.n_xpoints = 100;
            nbdopt.time = t;
            nbdopt.train.trial_idcs = 1;
            nbdopt.train.n_trials = numel(nbdopt.train.trial_idcs);
            nbdopt.test.trial_idcs = 2;
            nbdopt.test.n_trials = numel(nbdopt.test.trial_idcs);
            nbdopt.assumepoissonmdl = true;
            nbdopt.verbose = false;
            
            tic
            P_tR(:,:,ee,kk,bb) = naivebayestimedecoder(R,nbdopt);
            toc
            
            %% compute decoding performance statistics
            
            % iterate through time points
            for tt = 1 : T
                P_tR_mu(tt,ee,kk,bb) = P_tR(tt,:,ee,kk,bb) * t';
                P_tR_sd(tt,ee,kk,bb) = sqrt(...
                    P_tR(tt,:,ee,kk,bb) * (P_tR_mu(tt,ee,kk,bb) - t') .^ 2);
            end
        end
    end
end

%% choice of average & error functions
avgfun = @(x,d)nanmean(x,d);
errfun = @(x,d)nanstd(x,0,d);

%% plot cluster-split posterior averages
figure(...
    'name','cluster-split posterior averages',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(n_epochs,n_clusters);
for ee = 1 : n_epochs
    for kk = 1 : n_clusters
        sp_idx = kk + (ee - 1) * n_clusters;
        sps(ee,kk) = subplot(n_epochs,n_clusters,sp_idx);
        if ee == 1
            title(sps(ee,kk),cluster_labels{kk});
        end
    end
    xlabel(sps(ee,:),...
        sprintf('Time since %s onset (ms)',upper(epochs{ee})));
    ylabel(sps(ee,:),...
        sprintf('Decoded time since %s onset (ms)',upper(epochs{ee})));
end
set(sps,...
    axesopt.default,...
    'xlim',[ti,tf],...
    'xtick',unique([[ti,tf]';[ti,tf]';0;t_set]),...
    'ylim',[ti,tf],...
    'ytick',unique([[ti,tf]';[ti,tf]';0;t_set]));
linkaxes(sps);

clims = quantile(P_tR,[0,.999],'all')';

% iterate through epochs
for ee = 1 : n_epochs
    
    % iterate through clusters
    for kk = 1 : n_clusters
        
        p_tR = squeeze(avgfun(P_tR(:,:,ee,kk,:),5));
        imagesc(sps(ee,kk),[ti,tf],[ti,tf],p_tR',clims);
        plot(sps(ee,kk),xlim,ylim,'--w');
    end
end

%% plot superimposed cluster-split posterior averages

% figure initialization
fig = figure(figopt,...
    'position',[250,350,900,420],...
    'name','ramps_superimposed_posterior_averages_s1s2',...
    'numbertitle','off');

% axes initialization
n_sps = n_epochs;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(1,n_epochs,ii);
    title(sps(ii),epochs{ii});
    xlabel(sps(ii),sprintf('Time since %s onset (ms)',upper(epochs{ii})));
    ylabel(sps(ii),sprintf('Decoded time since %s onset (ms)',upper(epochs{ii})));
end
set(sps,...
    axesopt.default,...
    'xlim',[ti,tf],...
    'xtick',unique([[ti,tf]';[ti,tf]';0;t_set]),...
    'ylim',[ti,tf],...
    'ytick',unique([[ti,tf]';[ti,tf]';0;t_set]));

% iterate through epochs
for ee = 1 : n_epochs
    
    % convert from tensor to rgb
    P_tR_avg = squeeze(avgfun(P_tR(:,:,ee,:,:),5));
    P = tensor2rgb(permute(P_tR_avg,[2,1,3]),cluster_clrs);
    imagesc(sps(ee),[ti,tf],[ti,tf],P);
    
    % zero lines
    plot(sps(ee),xlim,[1,1]*0,':k');
    
    % plot identity line
    plot(sps(ee),[ti,t_set(end)],[ti,t_set(end)],':k');
    
    % iterate through cluster
    for kk = 1 : n_clusters
        pthat_avg = squeeze(avgfun(P_tR_mu(:,ee,kk,:),4));
        pthat_err = squeeze(errfun(P_tR_mu(:,ee,kk,:),4));
        plot(sps(ee),t,pthat_avg,...
            'color',cluster_clrs(kk,:),...
            'linewidth',1.5);
    end
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% decoding accuracy in ramps & non-ramps across clusters

% figure initialization
fig = figure(figopt,...
    'position',[100,50,440,420],...
    'name','ramps_decoding_accuracy_s1s2');

% axes initialization
xxtick = unique((1:n_epochs*2)+[-1;0;1]*.05*n_epochs*2);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:n_epochs)) = {''};
xxticklabel(ismember(xxtick,1:n_epochs)) = upper(epochs);
axes(axesopt.default,...
    'plotboxaspectratio',[2.5,1,1],...
    'color','none',...
    'xlim',[1,n_epochs*2]+[-1,1]*.1*n_epochs*2,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'clipping','off',...
    'layer','bottom');
xlabel('Epochs');
ylabel('Error (ms)');

% offset between ramps and non-ramps
xoffsets = [-1,1] * .15;

% choice of accuracy function
accuracyfun = @(x,d) nanmean(abs(x-t'),d);

% choice of average and error functions
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);

% preallocation
accuracy = nan(T,n_epochs,n_clusters);
% accuracy = nan(n_boots,n_epochs,n_clusters);

% iterate through epochs
for ee = 1 : n_epochs
    
    % iterate through clusters
    for kk = 1 : n_clusters
        accuracy(:,ee,kk) = accuracyfun(P_tR_mu(:,ee,kk,:),4);
        avg = avgfun(accuracy(:,ee,kk));
        err = errfun(accuracy(:,ee,kk));
        errorbar(ee+xoffsets(1+iseven(kk)),avg,...
            err(1),err(2),...
            'color','k',...
            'marker','o',...
            'markersize',7.5,...
            'markeredgecolor','k',...
            'markerfacecolor',cluster_clrs(kk,:),...
            'linewidth',1.5,...
            'capsize',0);
    end
end

% update axes
yymax = ceil(max(ylim)/10) * 10;
yylim = [0,yymax];
yytick = linspace(yylim(1),yylim(2),2);
yyticklabel = num2cell(yytick);
yyticklabel(~ismember(yytick,yylim)) = {''};
set(gca,...
    'ylim',yylim + [-1,1] * .05 * 2.5 * range(yylim),...
    'ytick',yytick,...
    'yticklabel',yyticklabel); % {'0',''});

% iterate through epochs
for ee = 1 : n_epochs
    xx = [-1,1] * .5 / 3 + ee;
    yy = [1,1] * max(yylim);
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    pval = kruskalwallis([accuracy(:,ee,1),accuracy(:,ee,2)],[],'off');
    if pval < .01
        test_str = '**';
    elseif pval < .05
        test_str = '*';
    else
        test_str = 'n.s.';
    end
    text(mean(xx),mean(yy)-.025/2.5*range(ylim),test_str,...
        'color','k',...
        'fontsize',16,...
        'horizontalalignment','center',...
        'verticalalignment','bottom');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% decoding precision in ramps & non-ramps across clusters

% figure initialization
fig = figure(figopt,...
    'position',[100,350,440,420],...
    'name','ramps_decoding_precision_s1s2');

% axes initialization
xxtick = unique((1:n_epochs*2)+[-1;0;1]*.05*n_epochs*2);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:n_epochs)) = {''};
xxticklabel(ismember(xxtick,1:n_epochs)) = upper(epochs);
axes(axesopt.default,...
    'plotboxaspectratio',[2.5,1,1],...
    'color','none',...
    'xlim',[1,n_epochs*2]+[-1,1]*.1*n_epochs*2,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'clipping','off',...
    'layer','bottom');
xlabel('Epochs');
ylabel('SD (ms)');

% offset between ramps and non-ramps
xoffsets = [-1,1] * .15;

% choice of average and error functions
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);

% preallocation
precision = nan(T,n_epochs,n_clusters);
% precision = nan(n_boots,n_epochs,n_clusters);

% iterate through epochs
for ee = 1 : n_epochs
    
    % iterate through clusters
    for kk = 1 : n_clusters
        precision(:,ee,kk) = nanmean(P_tR_sd(:,ee,kk,:),4);
        avg = avgfun(precision(:,ee,kk));
        err = errfun(precision(:,ee,kk));
        errorbar(ee+xoffsets(1+iseven(kk)),avg,...
            err(1),err(2),...
            'color','k',...
            'marker','o',...
            'markersize',7.5,...
            'markeredgecolor','k',...
            'markerfacecolor',cluster_clrs(kk,:),...
            'linewidth',1.5,...
            'capsize',0);
    end
end

% update axes
yymax = ceil(max(ylim)/10) * 10;
yylim = [0,yymax];
yytick = linspace(yylim(1),yylim(2),2);
yyticklabel = num2cell(yytick);
yyticklabel(~ismember(yytick,yylim)) = {''};
set(gca,...
    'ylim',yylim + [-1,1] * .05 * 2.5 * range(yylim),...
    'ytick',yytick,...
    'yticklabel',yyticklabel); % {'0',''});

% iterate through epochs
for ee = 1 : n_epochs
    xx = [-1,1] * .5 / 3 + ee;
    yy = [1,1] * max(yylim);
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    pval = kruskalwallis([precision(:,ee,1),precision(:,ee,2)],[],'off');
    if pval < .01
        test_str = '**';
    elseif pval < .05
        test_str = '*';
    else
        test_str = 'n.s.';
    end
    text(mean(xx),mean(yy)-.025/2.5*range(ylim),test_str,...
        'color','k',...
        'fontsize',16,...
        'horizontalalignment','center',...
        'verticalalignment','bottom');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end