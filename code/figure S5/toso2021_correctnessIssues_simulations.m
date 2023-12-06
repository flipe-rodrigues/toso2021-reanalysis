%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% simulation settings
N = 20;
T = n_tbins;
P = size(s_pairset,1);
S = size(s_pairset,2);
K = 100;

%% ROI & time settings
si_roi = [0,t_set(end)];
ti = si_roi(1);
tf = si_roi(2);
t = linspace(ti,tf*2,T);

%% model parameters
mus = [linspace(ti,tf,N/2),linspace(ti,tf,N/2)];
sigma = .05 * range(si_roi);
x = normpdf(t',mus,sigma);

%% temporal scale settings
tsi = 0;
tsf = 2;
ts_support = linspace(tsi,tsf,1e2);
ts_mu = 1;
ts_sigma = .05;
ts_pd = truncate(makedist('normal',...
    'mu',ts_mu,...
    'sigma',ts_sigma),0,inf);
ts_pdf = ts_pd.pdf(ts_support);
figure; plot(ts_support,ts_pdf);

%% color settings
slow_clr = [.0,.4,.95];
avg_clr = [1,1,1] * .0;
fast_clr = [.95,.25,.35];
clrmap_res = 1e3;
clrmap = colorlerp([slow_clr;avg_clr;fast_clr],clrmap_res);

%% generate fake data

% sample stimulus durations
t1_samples = datasample(t1(valid_flags),K);
t2_samples = datasample(t2(valid_flags),K);
si_samples = datasample([t1(valid_flags),t2(valid_flags)],K);

% sample temporal scaling factors
ts = random(ts_pd,K,S);

% preallocation
X = nan(T,N,K,S);

% iterate through neurons
for nn = 1 : N
         
    % iterate through trials
    for kk = 1 : K
        
        % iterate through stimuli
        for ss = 1 : S
            xs = interp1(t,x(:,nn),t * ts(kk,ss));
            si_flags = t <= si_samples(kk,ss);
            X(si_flags,nn,kk,ss) = xs(si_flags);
        end
    end
end

% normalization
X = (X - min(X,[],'all')) / range(X,'all');
x = (x - min(x,[],'all')) / range(x,'all');

%%

% figure initialization
fig = figure(figopt,...
    'name','correctness_issues');

% axes initialization
n_rows = 1;
n_cols = S;
n_sps = n_rows * n_cols;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(n_rows,n_cols,ii);
end
set(sps,axesopt.default,...
    'xlim',si_roi,...
    'ylim',[0,(N+1)*1.25],...
    'xtick',unique([si_roi,t_set']),...
    'ytick',unique([si_roi,t_set']),...
    'ycolor','none',...
    'clipping','off');
xlabel('Time (ms) X_i');

%
trial_clrs = colorlerp([slow_clr;avg_clr;fast_clr],K);

% iterate through trials
for kk = 1 : K
    plot(t,X(:,:,kk,1).*(1:N <= N/2)+(1:N)*1.05,'k');
    plot(t+T*3,X(:,:,kk,1).*(1:N > N/2)+(1:N)*1.05,'k');
end


% iterate through trials
for kk = 1 : K
    plot(t,X(:,:,kk,1).*(1:N <= N/2)+(1:N)*1.05,'k');
    plot(t+T*3,X(:,:,kk,1).*(1:N > N/2)+(1:N)*1.05,'k');
end
axis tight;

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%%

% figure initialization
fig = figure(figopt,...
    'name','correctness_issues');

% axes initialization
axes(axesopt.default,...
    'xlim',si_roi,...
    'ylim',[0,(N+1)*1.25],...
    'xtick',unique([si_roi,t_set']),...
    'ytick',unique([si_roi,t_set']),...
    'ycolor','none',...
    'clipping','off');
xlabel('Time (ms) X_i');

nah = sort(random(ts_pd,1,T));
for nn = 1 : N
    spacer = nn * 1.25;
    xpatch = [t,fliplr(t)];
    ypatch = [interp1(t,x(:,nn)',t*.75),fliplr(interp1(t,x(:,nn)',t*1.25))] + spacer;
    ypatch(isnan(ypatch)) = spacer;
    cpatch = [nah,fliplr(nah)];
    patch(xpatch,ypatch,ypatch-spacer,...
        'edgecolor','interp',...
        'facecolor','none',...
        'facealpha',1,...
        'linewidth',1.5,...
        'linestyle','-');
end

% save figure
if want2save
%     svg_file = fullfile(panel_path,[fig.Name,'.svg']);
%     print(fig,svg_file,'-dsvg','-painters');
end

%%
figure;
hold on;
nah = sort(random(ts_pd,1,T));
for nn = 1 : N
    x_spacer = (nn > N/2) * T * 3;
    y_spacer = nn * 1.25;
    
    y = interp1(t,x(:,nn),t*.75);
    xpatch = [fliplr(t),t] + x_spacer;
    ypatch = [zeros(1,T),y] + y_spacer;
    nan_flags = isnan(ypatch);
    patch(xpatch(~nan_flags),ypatch(~nan_flags),slow_clr,...
        'facealpha',1,...
        'edgecolor','none');
    plot(t + x_spacer,y + y_spacer,...
        'color','k',...
        'linewidth',1.5);

    y = interp1(t,x(:,nn),t*1.25);
    xpatch = [fliplr(t),t] + x_spacer;
    ypatch = [zeros(1,T),y] + y_spacer;
    nan_flags = isnan(ypatch);
    patch(xpatch(~nan_flags),ypatch(~nan_flags),fast_clr,...
        'facealpha',1,...
        'edgecolor','none');
    plot(t + x_spacer,y + y_spacer,...
        'color','k',...
        'linewidth',1.5);
    
    y = interp1(t,x(:,nn),t*1);
    xpatch = [fliplr(t),t] + x_spacer;
    ypatch = [zeros(1,T),y] + y_spacer;
    nan_flags = isnan(ypatch);
    patch(xpatch(~nan_flags),ypatch(~nan_flags),avg_clr,...
        'facealpha',1,...
        'edgecolor','none');
    plot(t + x_spacer,y + y_spacer,...
        'color','k',...
        'linewidth',1.5);
end
axis tight;


return;

%% compute cross-trial averages
x = nanmean(X,3);
r = nanmean(X,3);

%% compute correlation coefficients
rhos_monotonocity_onset = nan(N,1);
rhos_monotonocity_offset = nan(N,1);

% iterate through neurons
for nn = 1 : N
    rhos_monotonocity_onset(nn) = ...
        corr((1:T)',r(roi_onset_flags,nn));
    rhos_monotonocity_offset(nn) = ...
        corr((1:T)',r(roi_offset_flags,nn));
end

%% perform linear regression
mdls_onset = cell(N,1);
mdls_offset = cell(N,1);
pvals_monotonocity_onset = nan(N,1);
pvals_monotonocity_offset = nan(N,1);
betas_monotonocity_onset = nan(N,1);
betas_monotonocity_offset = nan(N,1);

% iterate through neurons
for nn = 1 : N
    mdls_onset{nn} = fitlm(1:T,r(roi_onset_flags,nn));
    mdls_offset{nn} = fitlm(1:T,r(roi_offset_flags,nn));
    p_onset = polyfit(1:T,r(roi_onset_flags,nn)',1);
    p_offset = polyfit(1:T,r(roi_offset_flags,nn)',1);
    pvals_monotonocity_onset(nn) = ...
        mdls_onset{nn}.Coefficients.pValue(end);
    pvals_monotonocity_offset(nn) = ...
        mdls_offset{nn}.Coefficients.pValue(end);
    betas_monotonocity_onset(nn) = p_onset(1);
    betas_monotonocity_offset(nn) = p_offset(1);
end

%% stereotypy assessment
k = floor(P / P);
shuffled_trial_idcs = randperm(P,P);

% preallocation
r_partitions_onset = nan(T,N,P);
r_partitions_offset = nan(T,N,P);
rhos_stereotypy_onset = nan(N,P);
rhos_stereotypy_offset = nan(N,P);
pvals_stereotypy_onset = nan(N,P);
pvals_stereotypy_offset = nan(N,P);

% iterate through partitions
for pp = 1 : P
    partition_idcs = shuffled_trial_idcs((1 : k) + (pp - 1) * k);
    r_partitions_onset(:,:,pp) = ...
        nanmean(X(roi_onset_flags,:,partition_idcs),3);
    r_partitions_offset(:,:,pp) = ...
        nanmean(X(roi_offset_flags,:,partition_idcs),3);
end

% compute reference
r_ref_onset = nanmean(r_partitions_onset,3);
r_ref_offset = nanmean(r_partitions_offset,3);

% iterate through partitions
for pp = 1 : P
    
    % iterate through neurons
    for nn = 1 : N
        [rhos_onset,pvals_onset] = ...
            corrcoef(r_ref_onset(:,nn),r_partitions_onset(:,nn,pp));
        [rhos_offset,pvals_offset] = ...
            corrcoef(r_ref_offset(:,nn),r_partitions_offset(:,nn,pp));
        rhos_stereotypy_onset(nn,pp) = rhos_onset(1,2);
        rhos_stereotypy_offset(nn,pp) = rhos_offset(1,2);
        pvals_stereotypy_onset(nn,pp) = pvals_onset(1,2);
        pvals_stereotypy_offset(nn,pp) = pvals_offset(1,2);
    end
end

%% neuron selection

% flag "monotonic" neurons
monotonicity_flags_onset = ...
    abs(rhos_monotonocity_onset) > rho_monotonocity_cutoff & ...
    abs(betas_monotonocity_onset) > beta_monotonocity_cutoff & ...
    pvals_monotonocity_onset <= pval_monotonocity_cutoff;
monotonicity_flags_offset = ...
    abs(rhos_monotonocity_offset) > rho_monotonocity_cutoff & ...
    abs(betas_monotonocity_offset) > beta_monotonocity_cutoff & ...
    pvals_monotonocity_offset <= pval_monotonocity_cutoff;

% flag stereotypical neurons
stereotypy_flags_onset = ...
    mean(rhos_stereotypy_onset > rho_stereotypy_cutoff,2) == 1 & ...
    mean(pvals_stereotypy_onset < pval_stereotypy_cutoff,2) == 1;
stereotypy_flags_offset = ...
    mean(rhos_stereotypy_offset > rho_stereotypy_cutoff,2) == 1 & ...
    mean(pvals_stereotypy_offset < pval_stereotypy_cutoff,2) == 1;

% flag "ramping" neurons
ramp_flags_onset = ...
    monotonicity_flags_onset & ...
    stereotypy_flags_onset;
ramp_flags_offset = ...
    monotonicity_flags_offset & ...
    stereotypy_flags_offset;
ramp_flags = ...
    ramp_flags_onset | ...
    ramp_flags_offset;

%% naive bayes decoder
train_flags = ismember(1:P,randperm(P,round(P/2)));
tensor_ramp = cat(3,...
    nanmean(X(roi_decoding_flags,ramp_flags,train_flags),3),...
    nanmean(X(roi_decoding_flags,ramp_flags,~train_flags),3));
tensor_non = cat(3,...
    nanmean(X(roi_decoding_flags,~ramp_flags,train_flags),3),...
    nanmean(X(roi_decoding_flags,~ramp_flags,~train_flags),3));

% enforce equal numbers of neurons on both clusters
n_ramp = sum(ramp_flags);
n_non = sum(~ramp_flags);
n = min(n_ramp,n_non);
all_ramp_idcs = all_idcs(ramp_flags);
all_non_idcs = all_idcs(~ramp_flags);
selected_ramp_idcs = sort(randperm(n_ramp,N_clus));
selected_non_idcs = sort(randperm(n_non,N_clus));
selected_idcs = sort([...
    all_ramp_idcs(selected_ramp_idcs),...
    all_non_idcs(selected_non_idcs)]);
tensor_ramp = tensor_ramp(:,selected_ramp_idcs,:);
tensor_non = tensor_non(:,selected_non_idcs,:);

% decoding options
opt = struct();
opt.n_xpoints = 100;
opt.time = t_roi_decoding;
opt.train.trial_idcs = 1;
opt.train.n_trials = numel(opt.train.trial_idcs);
opt.test.trial_idcs = 2;
opt.test.n_trials = numel(opt.test.trial_idcs);
opt.assumepoissonmdl = true;
opt.verbose = false;

% preallocation
P_tR_ramp = naivebayestimedecoder(tensor_ramp,opt);
P_tR_non = naivebayestimedecoder(tensor_non,opt);

%% compute decoding statistics

% preallocation
mu_ramp = nan(T,1);
mu_non = nan(T,1);
sd_ramp = nan(T,1);
sd_non = nan(T,1);

% iterate through time points
for tt = 1 : T
    mu_ramp(tt) = P_tR_ramp(tt,:) * t_roi_decoding';
    mu_non(tt) = P_tR_non(tt,:) * t_roi_decoding';
    sd_ramp(tt) = sqrt(P_tR_ramp(tt,:) * (mu_ramp(tt) - t_roi_decoding') .^ 2);
    sd_non(tt) = sqrt(P_tR_non(tt,:) * (mu_non(tt) - t_roi_decoding') .^ 2);
end

% compute posterior median
median_flags_ramp = [false(T,1),diff(cumsum(P_tR_ramp,2) > .5,1,2) == 1];
median_flags_non = [false(T,1),diff(cumsum(P_tR_non,2) > .5,1,2) == 1];
[~,median_idcs_ramp] = max(median_flags_ramp,[],2);
[~,median_idcs_non] = max(median_flags_non,[],2);
med_ramp = t_roi_decoding(median_idcs_ramp);
med_non = t_roi_decoding(median_idcs_non);

% compute posterior IQR
q25_flags_ramp = [false(T,1),diff(cumsum(P_tR_ramp,2) > .25,1,2) == 1];
q75_flags_ramp = [false(T,1),diff(cumsum(P_tR_ramp,2) > .75,1,2) == 1];
q25_flags_non = [false(T,1),diff(cumsum(P_tR_non,2) > .25,1,2) == 1];
q75_flags_non = [false(T,1),diff(cumsum(P_tR_non,2) > .75,1,2) == 1];
[~,q25_idcs_ramp] = max(q25_flags_ramp,[],2);
[~,q75_idcs_ramp] = max(q75_flags_ramp,[],2);
[~,q25_idcs_non] = max(q25_flags_non,[],2);
[~,q75_idcs_non] = max(q75_flags_non,[],2);
iqr_ramp = t_roi_decoding(q75_idcs_ramp) - t_roi_decoding(q25_idcs_ramp);
iqr_non = t_roi_decoding(q75_idcs_non) - t_roi_decoding(q25_idcs_non);

% compute MAP
[~,mode_idcs_ramp] = max(P_tR_ramp,[],2);
[~,mode_idcs_non] = max(P_tR_non,[],2);
map_ramp = t_roi_decoding(mode_idcs_ramp);
map_non = t_roi_decoding(mode_idcs_non);

%% store current simulation
P_TR_RAMP(:,:,ss,mm) = P_tR_ramp;
P_TR_NON(:,:,ss,mm) = P_tR_non;
MAP_ramp(:,ss,mm) = map_ramp;
MAP_non(:,ss,mm) = map_non;
MU_ramp(:,ss,mm) = mu_ramp;
MU_non(:,ss,mm) = mu_non;
MED_ramp(:,ss,mm) = med_ramp;
MED_non(:,ss,mm) = med_non;
SD_ramp(:,ss,mm) = sd_ramp;
SD_non(:,ss,mm) = sd_non;
IQR_ramp(:,ss,mm) = iqr_ramp;
IQR_non(:,ss,mm) = iqr_non;
MUS(:,ss,mm) = mus(selected_idcs);
SIGMAS(:,ss,mm) = sigma(selected_idcs);
GAMMAS(:,ss,mm) = gammas(selected_idcs);
ETAS(:,ss,mm) = etas(selected_idcs);
MEAN_FR(:,ss,mm) = mean(r(:,selected_idcs),1);
SELECTED_RAMP_FLAGS(:,ss,mm) = ramp_flags(selected_idcs);
ALL_RAMP_FLAGS(:,ss,mm) = ramp_flags;
P_RAMP(ss,mm) = nanmean(ramp_flags);

%% compute analysis metrics

% compute temporal tuning
tuning = t * ((r - nanmin(r)) ./ nansum(r - nanmin(r)));

% compute firing rate range
fr_range = range(r);

% preallocation
stereotypy = nan(N,1);

% iterate through neurons
for nn = 1 : N
    
    % compute stereotypy
    train_flags = ismember(1:P,randperm(P,P/2));
    rho = corrcoef(...
        nanmean(X(:,nn,train_flags),3),...
        nanmean(X(:,nn,~train_flags),3));
    stereotypy(nn) = rho(1,2);
end

% store metrics
TUNING(:,ss,mm) = tuning(selected_idcs);
FR_RANGE(:,ss,mm) = fr_range(selected_idcs);
STEREOTYPY(:,ss,mm) = stereotypy(selected_idcs);

%% model selection
model2plot = M;
ramps2plot_flags = SELECTED_RAMP_FLAGS(:,:,model2plot);

%% parameter selection

% preallocation
params2plot = struct();

% assignment
params2plot.mu = MUS(:,:,model2plot);
params2plot.gamma = GAMMAS(:,:,model2plot);
params2plot.eta = ETAS(:,:,model2plot);

% parse selected parameters
params2plot_labels = fieldnames(params2plot);
n_params2plot = numel(params2plot_labels);

%% metric selection

% preallocation
metrics2plot = struct();

% assignment
metrics2plot.tuning = TUNING(:,:,model2plot);
metrics2plot.frrange = FR_RANGE(:,:,model2plot);
metrics2plot.stereotypy = STEREOTYPY(:,:,model2plot);

% parse selected metrics
metrics2plot_labels = fieldnames(metrics2plot);
n_metrics2plot = numel(metrics2plot_labels);

%% parameter distributions

% figure initialization
fig = figure(figopt,...
    'position',[200,200,560,450],...
    'name','ramp_parameter_distributions');

% axes initialization
n_sps = n_params2plot;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(n_sps,1,ii);
    ylabel(sps(ii),'PDF');
end
set(sps,axesopt.default,...
    'plotboxaspectratio',[5,1,1],...
    'ticklength',axesopt.default.ticklength,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'ytick',0,...
    'clipping','off');
xlabel(sps(1),'\mu (ms)');
xlabel(sps(2),'\gamma');
xlabel(sps(3),'\eta (ms)');

% preallocation
bounds = struct();
edges = struct();

% bin settings
n_bins = 30;
bounds.mu = [ti,tf];
bounds.gamma = [0,25];
bounds.eta = [0,range(si_roi)*3/4];

% iterate through selected parameters
for ii = 1 : n_params2plot
    param = params2plot_labels{ii};
    
    % compute parameter distributions
    edges.(param) = linspace(bounds.(param)(1),bounds.(param)(2),n_bins);
    counts_all = histcounts(...
        params2plot.(param),edges.(param));
    counts_non = histcounts(...
        params2plot.(param)(~ramps2plot_flags),edges.(param));
    counts_ramp = histcounts(...
        params2plot.(param)(ramps2plot_flags),edges.(param));
    
    % plot distribution
    histogram(sps(ii),...
        'binedges',edges.(param),...
        'bincounts',counts_all,...
        'facecolor','w',...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    stairs(sps(ii),edges.(param),[counts_all,0],...
        'color','k',...
        'linewidth',1.5);
    histogram(sps(ii),...
        'binedges',edges.(param),...
        'bincounts',counts_non,...
        'facecolor',ramp_clrs(2,:),...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    stairs(sps(ii),edges.(param),[counts_non,0],...
        'color','k',...
        'linewidth',1.5);
    histogram(sps(ii),...
        'binedges',edges.(param),...
        'bincounts',counts_ramp,...
        'facecolor',ramp_clrs(1,:),...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    histogram(sps(ii),...
        'binedges',edges.(param),...
        'bincounts',counts_non,...
        'facecolor',ramp_clrs(2,:),...
        'edgecolor','none',...
        'facealpha',.5,...
        'linewidth',1.5);
    stairs(sps(ii),edges.(param),[counts_ramp,0],...
        'color','k',...
        'linewidth',1.5);
    
    % update axes
    set(sps(ii),...
        'xlim',bounds.(param),...
        'xtick',bounds.(param));
end

% update temporal tuning axes
set(sps(1),...
    'xlim',bounds.mu,...
    'xtick',unique([roi_onset,si_roi,roi_offset,t_set']));

% within-epoch unimodality assessment of temporal tuning
rois = [roi_onset; si_roi; roi_offset];
n_rois = size(rois,1);
for rr = 1 : n_rois
    roi_mus = params2plot.mu(ramps2plot_flags);
    mu_flags = ...
        roi_mus >= rois(rr,1) & ...
        roi_mus < rois(rr,2);
    [~,pval] = diptest(roi_mus(mu_flags));
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% metric distributions

% figure initialization
fig = figure(figopt,...
    'position',[200,200,560,450],...
    'name','ramp_metric_distributions');

% axes initialization
n_sps = n_metrics2plot;
sps = gobjects(n_sps,1);
for ii = 1 : n_sps
    sps(ii) = subplot(n_sps,1,ii);
    ylabel(sps(ii),'PDF');
end
set(sps,axesopt.default,...
    'plotboxaspectratio',[5,1,1],...
    'ticklength',axesopt.default.ticklength,...
    'xlimspec','tight',...
    'ylimspec','tight',...
    'ytick',0,...
    'clipping','off');
xlabel(sps(1),'Temporal tuning (ms)');
xlabel(sps(2),'Firing rate range (Hz)');
xlabel(sps(3),'Stereotypy coefficient');

% preallocation
bounds = struct();
edges = struct();

% bin settings
n_bins = 30;
bounds.tuning = si_roi;
bounds.frrange = [0,15];
bounds.stereotypy = [-1,1];

% iterate through selected metrics
for ii = 1 : n_metrics2plot
    metric = metrics2plot_labels{ii};
    
    % compute metric distributions
    edges.(metric) = linspace(bounds.(metric)(1),bounds.(metric)(2),n_bins);
    counts_all = histcounts(...
        metrics2plot.(metric),edges.(metric));
    counts_non = histcounts(...
        metrics2plot.(metric)(~ramps2plot_flags),edges.(metric));
    counts_ramp = histcounts(...
        metrics2plot.(metric)(ramps2plot_flags),edges.(metric));
    
    % plot distribution
    histogram(sps(ii),...
        'binedges',edges.(metric),...
        'bincounts',counts_all,...
        'facecolor','w',...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    stairs(sps(ii),edges.(metric),[counts_all,0],...
        'color','k',...
        'linewidth',1.5);
    histogram(sps(ii),...
        'binedges',edges.(metric),...
        'bincounts',counts_non,...
        'facecolor',ramp_clrs(2,:),...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    stairs(sps(ii),edges.(metric),[counts_non,0],...
        'color','k',...
        'linewidth',1.5);
    histogram(sps(ii),...
        'binedges',edges.(metric),...
        'bincounts',counts_ramp,...
        'facecolor',ramp_clrs(1,:),...
        'edgecolor','none',...
        'facealpha',1,...
        'linewidth',1.5);
    histogram(sps(ii),...
        'binedges',edges.(metric),...
        'bincounts',counts_non,...
        'facecolor',ramp_clrs(2,:),...
        'edgecolor','none',...
        'facealpha',.5,...
        'linewidth',1.5);
    stairs(sps(ii),edges.(metric),[counts_ramp,0],...
        'color','k',...
        'linewidth',1.5);
    
    % update axes
    set(sps(ii),...
        'xlim',bounds.(metric),...
        'xtick',bounds.(metric));
end

% update temporal tuning axes
set(sps(1),...
    'xlim',bounds.tuning,...
    'xtick',unique([roi_onset,si_roi,roi_offset,t_set']));

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% point estimate selection
pthat_str = 'MU';
pthat_ramp = eval([pthat_str,'_ramp']);
pthat_non = eval([pthat_str,'_non']);
pthat_all = eval([pthat_str,'_all']);
errhat_str = 'SD';
errhat_ramp = eval([errhat_str,'_ramp']);
errhat_non = eval([errhat_str,'_non']);
errhat_all = eval([errhat_str,'_all']);

%% decoding accuracy in ramps & non-ramps across models

% figure initialization
fig = figure(figopt,...
    'position',[100,50,440,420],...
    'name','ramps_decoding_accuracy');

% axes initialization
xxmax = M + 2;
xxtick = unique((1:xxmax)+[-1;0;1]*.05*xxmax);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:xxmax)) = {''};
xxticklabel(ismember(xxtick,1:M)) = model_labels;
xxticklabel(ismember(xxtick,M+[1,2])) = {'S1';'S2'};
axes(axesopt.default,...
    'plotboxaspectratio',[2.5,1,1],...
    'color','none',...
    'xlim',[1,xxmax]+[-1,1]*.1*xxmax,...
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
accuracyfun = @(x,d) nanmean(abs(x-t_roi_decoding'),d);

% choice of average and error functions
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);

% preallocation
ramp_sims = nan(S,M);
non_sims = nan(S,M);
% all_sims = nan(S,M);

% iterate through models
for mm = 1 : M
    ramp_sims(:,mm) = accuracyfun(pthat_ramp(:,:,mm),1);
    non_sims(:,mm) = accuracyfun(pthat_non(:,:,mm),1);
    ramp_avg = avgfun(ramp_sims(:,mm));
    non_avg = avgfun(non_sims(:,mm));
    ramp_err = errfun(ramp_sims(:,mm));
    non_err = errfun(non_sims(:,mm));
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
    'yticklabel',yyticklabel);

% iterate through models
for mm = 1 : M
    xx = [-1,1] * .5 / 3 + mm;
    yy = [1,1] * max(yylim);
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    [~,pval] = kstest2(ramp_sims(:,mm),non_sims(:,mm));
    pval = pval * M;
    if pval < .01
        test_str = '**';
        font_size = 16;
    elseif pval < .05
        test_str = '*';
        font_size = 16;
    else
        test_str = 'n.s.';
        font_size = axesopt.default.fontsize;
    end
    text(mean(xx),mean(yy)-.025/2.5*range(ylim),test_str,...
        'color','k',...
        'fontsize',font_size,...
        'horizontalalignment','center',...
        'verticalalignment','bottom');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% decoding precision in ramps & non-ramps across models

% figure initialization
fig = figure(figopt,...
    'position',[100,350,440,420],...
    'name','ramps_decoding_precision');

% axes initialization
xxmax = M + 2;
xxtick = unique((1:xxmax)+[-1;0;1]*.05*xxmax);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:xxmax)) = {''};
xxticklabel(ismember(xxtick,1:M)) = model_labels;
xxticklabel(ismember(xxtick,M+[1,2])) = {'S1';'S2'};
axes(axesopt.default,...
    'plotboxaspectratio',[2.5,1,1],...
    'color','none',...
    'xlim',[1,xxmax]+[-1,1]*.1*xxmax,...
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
ramp_sims = nan(S,M);
non_sims = nan(S,M);
% all_sims = nan(S,M);

% iterate through models
for mm = 1 : M
    ramp_sims(:,mm) = nanmean(errhat_ramp(:,:,mm),1);
    non_sims(:,mm) = nanmean(errhat_non(:,:,mm),1);
    ramp_avg = avgfun(ramp_sims(:,mm));
    non_avg = avgfun(non_sims(:,mm));
    ramp_err = errfun(ramp_sims(:,mm));
    non_err = errfun(non_sims(:,mm));
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
    'yticklabel',yyticklabel);

% iterate through models
for mm = 1 : M
    xx = [-1,1] * .5 / 3 + mm;
    yy = [1,1] * max(yylim);
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    [~,pval] = kstest2(ramp_sims(:,mm),non_sims(:,mm));
    pval = pval * M;
    if pval < .01
        test_str = '**';
        font_size = 16;
    elseif pval < .05
        test_str = '*';
        font_size = 16;
    else
        test_str = 'n.s.';
        font_size = axesopt.default.fontsize;
    end
    text(mean(xx),mean(yy)-.025/2.5*range(ylim),test_str,...
        'color','k',...
        'fontsize',font_size,...
        'horizontalalignment','center',...
        'verticalalignment','bottom');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% compare decoding from ramps & non-ramps side by side

% iterate through models
for mm = 1 : M
    figure(figopt,...
        'name',sprintf('posterior_means (model %i)',mm),...
        'numbertitle','off',...
        'windowstyle','docked');
    n_rows = 1;
    n_cols = 2;
    sps = gobjects(n_rows,n_cols);
    for rr = 1 : n_rows
        for cc = 1 : n_cols
            sp_idx = cc + (rr - 1) * n_cols;
            sps(rr,cc) = subplot(n_rows,n_cols,sp_idx);
            xlabel(sps(rr,cc),'Time (ms)');
            ylabel(sps(rr,cc),'Decoded time (ms)');
        end
    end
    set(sps,...
        axesopt.default,...
        'xlim',si_roi,...
        'ylim',si_roi,...
        'xdir','normal',...
        'ydir','normal',...
        'nextplot','add',...
        'plotboxaspectratio',[1,1,1]);
    linkaxes(sps);
    
    % ramping posteriors
    title(sps(1),sprintf('Ramping neurons (%.0f%%)',nanmean(P_RAMP(:,mm))*100));
    errorbar(sps(1),...
        t_roi_decoding,nanmean(pthat_ramp(:,:,mm),2),...
        nanmean(errhat_ramp(:,:,mm),2),...
        'color',ramp_clrs(1,:),...
        'linewidth',.1,...
        'capsize',0);
    plot(sps(1),t_roi_decoding,nanmean(pthat_ramp(:,:,mm),2),...
        'color','w',...
        'linewidth',1.5);
    
    % non-ramping posteriors
    title(sps(2),sprintf('Non-ramping neurons (%.0f%%)',(1-nanmean(P_RAMP(:,mm)))*100));
    errorbar(sps(2),...
        t_roi_decoding,nanmean(pthat_non(:,:,mm),2),...
        nanmean(errhat_non(:,:,mm),2),...
        'color',ramp_clrs(2,:),...
        'linewidth',.1,...
        'capsize',0);
    plot(sps(2),t_roi_decoding,nanmean(pthat_non(:,:,mm),2),...
        'color','w',...
        'linewidth',1.5);
    
    % annotate model parameters
    text(sps(1),.05,.95,model_labels{mm},...
        'horizontalalignment','left',...
        'verticalalignment','top',...
        'units','normalized');
end

%% compare decoding from ramps & non-ramps (superimposed)

% figure initialization
fig = figure(figopt,...
    'name',sprintf('superimposed_posterior_means (model %i)',M));

% axes initialization
axes(axesopt.default,...
    'xlim',si_roi,...
    'ylim',si_roi,...
    'xtick',unique([si_roi,t_set']),...
    'ytick',unique([si_roi,t_set']),...
    'clipping','off');
xlabel('Time (ms) X_i');
ylabel('Decoded time (ms) X_i');

% choice of average function
avgfun = @(x,d) nanmean(x,d);

% non-ramping posteriors
non_avg = avgfun(pthat_non(:,:,M),2);
non_err = [1,1] .* avgfun(errhat_non(:,:,M),2);
xpatch = [t_roi_decoding,fliplr(t_roi_decoding)];
ypatch = [non_avg-non_err(:,1);flipud(non_avg+non_err(:,2))];
patch(xpatch,ypatch,ramp_clrs(2,:),...
    'facealpha',1,...
    'edgecolor','none');

% ramping posteriors
ramp_avg = avgfun(pthat_ramp(:,:,M),2);
ramp_err = [1,1] .* avgfun(errhat_ramp(:,:,M),2);
xpatch = [t_roi_decoding,fliplr(t_roi_decoding)];
ypatch = [ramp_avg-ramp_err(:,1);flipud(ramp_avg+ramp_err(:,2))];
patch(xpatch,ypatch,ramp_clrs(1,:),...
    'facealpha',1,...
    'edgecolor','none');

% plot averages
plot(t_roi_decoding,non_avg,...
    'color',([1,1,1]+ramp_clrs(2,:))/2,...
    'linewidth',1.5);
plot(t_roi_decoding,ramp_avg,...
    'color',([1,1,1]+ramp_clrs(1,:))/2,...
    'linewidth',1.5);

% inset with example run
axes(axesopt.default,...
    axesopt.inset.nw,...
    'xaxislocation','bottom',...
    'xlim',si_roi,...
    'ylim',si_roi,...
    'xtick',si_roi,...
    'ytick',si_roi);

% draw example
eg_idx = randi(S);
eg_clr = [1,1,1] * .75;

% example non-ramping posteriors
non_avg_eg = pthat_non(:,eg_idx,M);
non_err_eg = [1,1] .* errhat_non(:,eg_idx,M);
xpatch = [t_roi_decoding,fliplr(t_roi_decoding)];
ypatch = [non_avg_eg-non_err_eg(:,1);flipud(non_avg_eg+non_err_eg(:,2))];
patch(xpatch,ypatch,eg_clr,...
    'facealpha',1,...
    'edgecolor','none');

% plot example average
plot(t_roi_decoding,non_avg_eg,...
    'color',([1,1,1]+eg_clr)/2,...
    'linewidth',1.5);

% plot reference line
plot(t_roi_decoding,t_roi_decoding,'--k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% <firing rate> distributions across models
figure(...
    'position',[1.0258e+03 41.8000 1.0224e+03 1.0288e+03]);
fr_edges = linspace(0,30,50);
for mm = 1 : M
    subplot(M,1,mm);
    set(gca,axesopt.default,...
        'plotboxaspectratio',[M,1,1]);
    xlabel('<Firing rate> (Hz)');
    ylabel('Count');
    mm_ramp_flags = SELECTED_RAMP_FLAGS(:,:,mm);
    mm_mean_fr = MEAN_FR(:,:,mm);
    histogram(mm_mean_fr(mm_ramp_flags),fr_edges,...
        'facecolor',ramp_clrs(1,:),...
        'linewidth',1.5);
    histogram(mm_mean_fr(~mm_ramp_flags),fr_edges,...
        'facecolor',ramp_clrs(2,:),...
        'linewidth',1.5);
    axis tight;
    yymax = max(ylim);
    plot([1,1]*mean(mm_mean_fr(mm_ramp_flags),[1,2]),[0,1]*yymax*1.05,...
        'color',ramp_clrs(1,:),...
        'linewidth',3);
    plot([1,1]*mean(mm_mean_fr(~mm_ramp_flags),[1,2]),[0,1]*yymax*1.05,...
        'color',ramp_clrs(2,:),...
        'linewidth',3);
    text(1,1,sprintf('<P(ramp)> = %.2f',mean(P_RAMP(:,mm))),...
        'horizontalalignment','right',...
        'verticalalignment','top',...
        'units','normalized');
end