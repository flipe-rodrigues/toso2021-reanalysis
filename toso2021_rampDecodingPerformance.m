%% initialization
if ~exist('data','var')
    toso2021_wrapper;
    close all;
end

%% bootstrap settings
n_boots = 25;

%% concatenation settings
n_concats_max = 2^7;

%% temporal smoothing kernel
gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);

%% time settings
dt = 75;
ti = 0;
tf = t_set(end);
t = ti : tf;
t_padded = ti - dt/2 : tf + dt/2;
T = numel(t);
t_idcs = 1 : T;
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

%% construct spike rate tensor (time X neurons X concatenations)

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
            R = R_all(:,idcs(cluster_flags),:);
            
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
sps = gobjects(n_clusters,1);
for ii = 1 : n_clusters
    sps(ii) = subplot(1,n_clusters,ii);
    title(sps(ii),cluster_labels{ii});
    xlabel(sps(ii),'Time since S_2 onset (ms)');
    ylabel(sps(ii),'Decoded time since S_2 onset (ms)');
end
set(sps,...
    axesopt.default,...
    'xlim',roi,...
    'xtick',unique([roi';roi';0;t_set]),...
    'ylim',roi,...
    'ytick',unique([roi';roi';0;t_set]));
linkaxes(sps);

clims = quantile(P_tR,[0,.999],'all')';

% iterate through clusters
for ii = 1 : n_clusters
    p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
    imagesc(sps(ii),roi,roi,p_cond',clims);
    plot(sps(ii),xlim,ylim,'--w');
end

%% plot superimposed cluster-split posterior averages
fig = figure(...
    figopt,...
    'name','superimposed_posterior_averages_rapms',...
    'numbertitle','off');
axes(...
    axesopt.default,...
    'xlim',roi,...
    'xtick',unique([roi';roi';0;t_set]),...
    'ylim',roi,...
    'ytick',unique([roi';roi';0;t_set]),...
    'xticklabelrotation',0,...
    'yticklabelrotation',0);
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% convert from tensor to rgb
P_tR_avg = squeeze(avgfun(P_tR,4));
P = tensor2rgb(permute(P_tR_avg,[2,1,3]),cluster_clrs);
imagesc(roi,roi,P);

% zero lines
plot(xlim,[1,1]*0,':k');

% plot identity line
plot([ti,t_set(end)],[ti,t_set(end)],':k');

% inset with pseudo colorbar
axes(...
    axesopt.default,...
    'position',[0.625,0.65,0.2583,0.2717],...
    'yaxislocation','right',...
    'xcolor','none',...
    'xlim',[0,1],...
    'ylim',[0,1],...
    'ytick',0,...
    'colormap',colorlerp(...
    [cluster_clrs(1,:);[1,1,1];cluster_clrs(end,:)],2^8));
ylabel('P(t|R)',...
    'verticalalignment','middle',...
    'rotation',-90);

% colorbar settings
clrbar_width = .05;

% iterate through clusters
for ii = 1 : n_clusters
    
    % patch pseudo-colorbar
    xpatch = (1 - clrbar_width * n_clusters) + ...
        clrbar_width * ((ii - 1) + [0,1,1,0]);
    ypatch = [0,0,1,1];
    patch(xpatch,ypatch,cluster_clrs(ii,:),...
        'edgecolor','none',...
        'linewidth',1.5);
end

% inset with cluster-split MAPs
axes(...
    axesopt.inset.se,...
    axesopt.default,...
    'box','off',...
    'color','w',...
    'xlim',[0,tf],...
    'xtick',unique([roi';roi';0;t_set]),...
    'ylim',[0,tf],...
    'ytick',unique([roi';roi';0;t_set]));

% iterate through cluster
for ii = 1 : n_clusters
    map_avg = squeeze(avgfun(P_tR_mu(:,ii,:),3));
    map_err = squeeze(errfun(P_tR_mu(:,ii,:),3));
    errorpatch(t,map_avg,map_err,cluster_clrs(ii,:),...
        'facealpha',.25);
    plot(t,map_avg,...
        'color',cluster_clrs(ii,:),...
        'linewidth',1.5);
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
    'name','ramps_decoding_accuracy');

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

% 
mm = 1;
ramp_sims = accuracyfun(P_tR_mu(:,1,:),3);
non_sims = accuracyfun(P_tR_mu(:,2,:),3);
ramp_avg = avgfun(ramp_sims);
non_avg = avgfun(non_sims);
ramp_err = errfun(ramp_sims);
non_err = errfun(non_sims);
plot(mm+xoffsets,[ramp_avg,non_avg],...
    'color','k',...
    'linewidth',1.5);
errorbar(mm+xoffsets(1),ramp_avg,...
    ramp_err(1),ramp_err(2),...
    'color','k',...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor',ramp_clrs(1,:),...
    'linewidth',1.5,...
    'capsize',0);
errorbar(mm+xoffsets(2),non_avg,...
    non_err(1),non_err(2),...
    'color','k',...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor',ramp_clrs(2,:),...
    'linewidth',1.5,...
    'capsize',0);

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

% stats
xx = [-1,1] * .5 / 3 + mm;
yy = [1,1] * max(yylim);
plot(xx,yy,...
    'color','k',...
    'linewidth',1.5,...
    'handlevisibility','off');
[~,pval] = kstest2(ramp_sims,non_sims);
pval = kruskalwallis([ramp_sims,non_sims],[],'off');
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

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% decoding precision in ramps & non-ramps across clusters

% figure initialization
fig = figure(figopt,...
    'position',[550,350,440,420],...
    'name','ramps_decoding_precision');

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

%
mm = 1;
ramp_sims = squeeze(nanmean(P_tR_sd(:,1,:),3));
non_sims = squeeze(nanmean(P_tR_sd(:,2,:),3));
ramp_avg = avgfun(ramp_sims);
non_avg = avgfun(non_sims);
ramp_err = errfun(ramp_sims);
non_err = errfun(non_sims);
plot(mm+xoffsets,[ramp_avg,non_avg],...
    'color','k',...
    'linewidth',1.5);
errorbar(mm+xoffsets(1),ramp_avg,...
    ramp_err(1),ramp_err(2),...
    'color','k',...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor',ramp_clrs(1,:),...
    'linewidth',1.5,...
    'capsize',0);
errorbar(mm+xoffsets(2),non_avg,...
    non_err(1),non_err(2),...
    'color','k',...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor',ramp_clrs(2,:),...
    'linewidth',1.5,...
    'capsize',0);

% update axes
yymax = ceil(max(ylim)/10) * 10;
yylim = [0,yymax];
yytick = linspace(yylim(1),yylim(2),2);
yyticklabel = num2cell(yytick);
yyticklabel(~ismember(yytick,yylim)) = {''};
set(gca,...
    'ylim',yylim + [-1,1] * .05 * 2.5 * range(yylim),...
    'ytick',yytick,...
    'yticklabel',yyticklabel);

%
xx = [-1,1] * .5 / 3 + mm;
yy = [1,1] * max(yylim);
plot(xx,yy,...
    'color','k',...
    'linewidth',1.5,...
    'handlevisibility','off');
[~,pval] = kstest2(ramp_sims,non_sims);
pval = kruskalwallis([ramp_sims,non_sims],[],'off');
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

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end