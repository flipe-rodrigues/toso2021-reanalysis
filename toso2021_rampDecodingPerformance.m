%% initialization
if ~exist('data','var')
    toso2021_wrapper;
    close all;
end

%% bootstrap settings
n_boots = 25;

%% temporal smoothing kernel
gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);

%% time settings
dt = 75;
ti = 0;
tf = t_set(end);
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
cluster_labels = {'ramp','non'};
cluster_idcs = cell2table({ramp_idcs.s2,nonramp_idcs.s2},...
    'variablenames',cluster_labels);
n_clusters = numel(cluster_labels);
min_cluster_size = inf;
for kk = 1 : n_clusters
    cluster = cluster_labels{kk};
    idcs = cluster_idcs.(cluster){:};
    cluster_idcs.(cluster) = {idcs(ismember(idcs,flagged_neurons))};
    min_cluster_size = min(min_cluster_size,numel(cluster_idcs.(cluster){:}));
end
cluster_clrs = ramp_clrs;

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
                R_all(:,nn,ee,tt+1) = r_mu;
            end
        end
    end
    
    %% decode time from each epoch & cluster
    
    % iterate through epochs
    for ee = 1 : n_epochs
        
        % iterate through clusters
        for kk = 1 : n_clusters
            cluster = cluster_labels{kk};
            idcs = cluster_idcs.(cluster){:};
            cluster_size = numel(idcs);
            cluster_flags = ismember(1:cluster_size,...
                randperm(cluster_size,min_cluster_size));
            R = R_all(:,idcs(cluster_flags),ee,:);
            
            %% naive bayes decoder
            nbdopt = struct();
            nbdopt.n_xpoints = 100;
            nbdopt.time = t;
            nbdopt.train.trial_idcs = 1;
            nbdopt.train.n_trials = numel(nbdopt.train.trial_idcs);
            nbdopt.test.trial_idcs = 2;
            nbdopt.test.n_trials = numel(nbdopt.test.trial_idcs);
            nbdopt.assumepoissonmdl = true;
            nbdopt.verbose = true;
            
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
n_sps = n_clusters * n_epochs;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(n_clusters,n_epochs,ii);
    xlabel(sps(ii),'Time since S_2 onset (ms)');
    ylabel(sps(ii),'Decoded time since S_2 onset (ms)');
end
set(sps,...
    axesopt.default,...
    'xlim',[ti,tf],...
    'xtick',unique([[ti,tf]';[ti,tf]';0;t_set]),...
    'ylim',[ti,tf],...
    'ytick',unique([[ti,tf]';[ti,tf]';0;t_set]));
title(sps(1),epochs{1});
title(sps(2),epochs{2});
linkaxes(sps);

clims = quantile(P_tR,[0,.999],'all')';

% iterate through clusters
for kk = 1 : n_clusters
    
    % iterate through epochs
    for ee = 1 : n_clusters
        sp_idx = ee + (kk - 1) * n_epochs;
        p_tR = squeeze(avgfun(P_tR(:,:,ee,kk,:),5));
        imagesc(sps(sp_idx),[ti,tf],[ti,tf],p_tR',clims);
        plot(sps(sp_idx),xlim,ylim,'--w');
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
    'ytick',unique([[ti,tf]';[ti,tf]';0;t_set]),...
    'xticklabelrotation',0,...
    'yticklabelrotation',0);

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
%         errorpatch(t,pthat_avg,pthat_err,cluster_clrs(kk,:),...
%             'facealpha',.25,...
%             'parent',sps(ee));
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
xxtick = unique((1:M)+[-1;0;1]*.05*M);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:M)) = {''};
xxticklabel(ismember(xxtick,1:M)) = model_labels;
axes(axesopt.default,...
    'plotboxaspectratio',[2.5,1,1],...
    'color','none',...
    'xlim',[1,M]+[-1,1]*.1*M,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'clipping','off',...
    'layer','bottom');
xlabel('Parameter range');
ylabel('Error (ms)');

% offset between ramps and non-ramps
xoffsets = [-1,1] * .15;

% choice of accuracy function
accuracyfun = @(x,d) nanmean(abs(x-t'),d);

% choice of average and error functions
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);

% preallocation
ramp_sims = nan(T,n_epochs);
non_sims = nan(T,n_epochs);

% iterate through epochs
for ee = 1 : n_epochs
    ramp_sims(:,ee) = accuracyfun(P_tR_mu(:,ee,1,:),4);
    non_sims(:,ee) = accuracyfun(P_tR_mu(:,ee,2,:),4);
    ramp_avg = avgfun(ramp_sims(:,ee));
    non_avg = avgfun(non_sims(:,ee));
    ramp_err = errfun(ramp_sims(:,ee));
    non_err = errfun(non_sims(:,ee));
    plot(ee+xoffsets,[ramp_avg,non_avg],...
        'color','k',...
        'linewidth',1.5);
    errorbar(ee+xoffsets(1),ramp_avg,...
        ramp_err(1),ramp_err(2),...
        'color','k',...
        'marker','o',...
        'markersize',7.5,...
        'markeredgecolor','k',...
        'markerfacecolor',ramp_clrs(1,:),...
        'linewidth',1.5,...
        'capsize',0);
    errorbar(ee+xoffsets(2),non_avg,...
        non_err(1),non_err(2),...
        'color','k',...
        'marker','o',...
        'markersize',7.5,...
        'markeredgecolor','k',...
        'markerfacecolor',ramp_clrs(2,:),...
        'linewidth',1.5,...
        'capsize',0);
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
    [~,pval] = kstest2(ramp_sims(:,ee),non_sims(:,ee));
    pval = kruskalwallis([ramp_sims(:,ee),non_sims(:,ee)],[],'off');
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
xxtick = unique((1:M)+[-1;0;1]*.05*M);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:M)) = {''};
xxticklabel(ismember(xxtick,1:M)) = model_labels;
axes(axesopt.default,...
    'plotboxaspectratio',[2.5,1,1],...
    'color','none',...
    'xlim',[1,M]+[-1,1]*.1*M,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'clipping','off',...
    'layer','bottom');
xlabel('Parameter range');
ylabel('SD (ms)');

% offset between ramps and non-ramps
xoffsets = [-1,1] * .15;

% choice of average and error functions
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);


% preallocation
ramp_sims = nan(T,n_epochs);
non_sims = nan(T,n_epochs);

% iterate through epochs
for ee = 1 : n_epochs
    ramp_sims(:,ee) = nanmean(P_tR_sd(:,ee,1,:),4);
    non_sims(:,ee) = nanmean(P_tR_sd(:,ee,2,:),4);
    ramp_avg = avgfun(ramp_sims(:,ee));
    non_avg = avgfun(non_sims(:,ee));
    ramp_err = errfun(ramp_sims(:,ee));
    non_err = errfun(non_sims(:,ee));
    plot(ee+xoffsets,[ramp_avg,non_avg],...
        'color','k',...
        'linewidth',1.5);
    errorbar(ee+xoffsets(1),ramp_avg,...
        ramp_err(1),ramp_err(2),...
        'color','k',...
        'marker','o',...
        'markersize',7.5,...
        'markeredgecolor','k',...
        'markerfacecolor',ramp_clrs(1,:),...
        'linewidth',1.5,...
        'capsize',0);
    errorbar(ee+xoffsets(2),non_avg,...
        non_err(1),non_err(2),...
        'color','k',...
        'marker','o',...
        'markersize',7.5,...
        'markeredgecolor','k',...
        'markerfacecolor',ramp_clrs(2,:),...
        'linewidth',1.5,...
        'capsize',0);
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
    [~,pval] = kstest2(ramp_sims(:,ee),non_sims(:,ee));
    pval = kruskalwallis([ramp_sims(:,ee),non_sims(:,ee)],[],'off');
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