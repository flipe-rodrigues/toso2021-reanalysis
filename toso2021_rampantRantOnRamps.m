%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% notes
% beta depends on time units (s or ms), and firing rate;
% rho and beta are redundant;
% stereoptypy is a circular criterion for the point of decodability

%% TODO:
% - stereotypy
% - model comparison with polynomial regression?

%% seed fixing
% rng(0);

%% tensor settings
T = 100;    % time points
N = 50;     % neurons
K = 100;    % trials
B = 100;    % bootstrap iterations

%% time settings
ti = 0;
tf = 1000;
t = ti : tf;
T = numel(t);
% t = linspace(ti,tf,T);
t_idcs = 1 : T;
t_units = 1e3;
dt = 75;
% dt = diff(t(1:2));

%% "ramping" criteria

% "monotonocity" criteria
rho_monotonocity_cutoff = .5;
beta_monotonocity_cutoff = .004; % .05;
pval_monotonocity_cutoff = .05;

% stereotypy criteria
rho_stereotypy_cutoff = .5;
pval_stereotypy_cutoff = .01;

%% smoothing kernel settings
kernel_peak_time = 50;
kernel = gammakernel('peakx',kernel_peak_time,'binwidth',dt);
t_padded = ti + kernel.paddx(1) : dt : tf + kernel.paddx(end);

%% color settings
clrs = cool(N);

%% example neurons with opposite gains
K_eg = 10;

% preallocation
X_eg = nan(T,2,K_eg);
R_eg = nan(T,2,K_eg);

% parameter choices
gammas_eg = [1,3] * 3;
% gammas_eg = [1,1] * mean(gammas_eg);
mus_eg = [1,1] * .1 * tf;
lambdas_eg = [0,0];
% lambdas_eg = [0,.15*tf];
sigmas_eg = [1,1] * .15 * (tf - ti);

% iterate through neurons
for nn = 1 : 2
    for kk = 1 : K_eg
        X_eg(:,nn,kk) = generativerate(...
            t,gammas_eg(nn),mus_eg(nn),lambdas_eg(nn),sigmas_eg(nn));
%         x_padded = generativerate(...
%             t_padded,gammas_eg(nn),mus_eg(nn),lambdas_eg(nn),sigmas_eg(nn));
%         dur = tf - ti - kernel.paddx(1) + kernel.paddx(2);
%         [n,ts] = poissonprocess(x_padded,dur/t_units);
%         spk_times = ts * t_units + ti + kernel.paddx(1);
%         spk_counts = histcounts(spk_times,'binedges',t_padded);
%         spk_counts = spk_counts / (dt/t_units);
%         R_eg(:,nn,kk) = conv(spk_counts,kernel.pdf,'valid');
%         R_eg(:,nn,kk) = spk_counts / (dt/t_units);        
        [n,ts] = poissonprocess(X_eg(:,nn,kk),(tf-ti)/t_units);
        spk_times = ts * t_units + ti;
        spk_counts = histcounts(spk_times,'binedges',[t,t(end)+1]);
        R_eg(:,nn,kk) = movsum(spk_counts,dt) / (dt/t_units);
    end
end

% compute cross-trial mean
r_eg = nanmean(R_eg,3);

% perform linear regression
pvals_monotonocity_eg = nan(2,1);
betas_monotonocity_eg = nan(2,1);

% iterate through neurons
for nn = 1 : 2
    mdl = fitlm(t_idcs,r_eg(:,nn));
    p = polyfit(t_idcs,r_eg(:,nn)',1);
    pvals_monotonocity_eg(nn) = mdl.Coefficients.pValue(end);
    betas_monotonocity_eg(nn) = p(1);
end

[~,idx] = max(r_eg(:,1));
figure(figopt);
axes(axesopt.default,...
    'xcolor','k',...
    'xtick',unique([ti,t(idx),tf]),...
    'xticklabel',{num2str(ti),'\mu_0',num2str(tf)},...
    'ylim',[0,5]+[-1,1]*.15*range([0,5]),...
    'ytick',[0,5],...
    'clipping','off',...
    'plotboxaspectratio',[3,1,1]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');

plot(t,r_eg(:,1),...
    'color','k',...
    'linewidth',1.5);
plot(t,r_eg(:,2),...
    'color','r',...
    'linewidth',1.5);

% annotate regression statistics
text(.5,1,...
    '$\gamma\sigma\lambda\phi\mu\sim N(\mu_0,\lambda)$',...
    'color','k',...
    'fontsize',12,...
    'interpreter','latex',...
    'horizontalalignment','center',...
    'units','normalized');
text(.95,.85,...
    sprintf('\\beta = %.2f, p-value = %.2f',...
    betas_monotonocity_eg(1),pvals_monotonocity_eg(1)),...
    'color','k',...
    'horizontalalignment','right',...
    'units','normalized');
text(.95,.95,...
    sprintf('\\beta = %.2f, p-value = %.2f',...
    betas_monotonocity_eg(2),pvals_monotonocity_eg(2)),...
    'color','r',...
    'horizontalalignment','right',...
    'units','normalized');

yspacer = 5;
figure(figopt);
axes(axesopt.default,...
    'xlim',[ti,tf],...
    ...'ylim',[1,K_eg*yspacer]+[-1,1]*.05*K_eg,...
    'ycolor','none',...
    'ylimspec','tight',...
    'clipping','off',...
    'plotboxaspectratio',[1,1,1]);
xlabel('Time (ms)');
ylabel('Trial #');

for kk = 1 : K_eg
    plot(t,R_eg(:,1,kk)+kk*yspacer,...
        'color','k',...
        'linewidth',1.5);
    plot(t,R_eg(:,2,kk)+kk*yspacer,...
        'color','r',...
        'linewidth',1.5);
end

ylim(ylim+[-1,1]*.05*range(ylim));

text(1,.95,...
    sprintf('\\gamma = %.2f, \\lambda = %.2f',...
    gammas_eg(1),lambdas_eg(1)),...
    'color','k',...
    'horizontalalignment','right',...
    'units','normalized');
text(1,.9,...
    sprintf('\\gamma = %.2f, \\lambda = %.2f',...
    gammas_eg(2),lambdas_eg(2)),...
    'color','r',...
    'horizontalalignment','right',...
    'units','normalized');
% imagesc(t,[1,K],nanmean(R_eg(:,1,:),3)')

%% bootstrapping

% preallocation
MU_ramp = nan(T,B);
MU_non = nan(T,B);
SD_ramp = nan(B,1);
SD_non = nan(B,1);
p_ramp = nan(B,1);
RAMP_FLAGS = nan(N,B);
GAMMAS = nan(N,B);
LAMBDAS = nan(N,B);

% iterate through bootstrap iterations
for bb = 1 : B
    progressreport(bb,B,'bootstrapping')
    
    %% model parameters
    gammas = exprnd(15,N,1);
    lambdas = exprnd(.2*(tf-ti),N,1);
    mus = linspace(ti,tf,N)';
    sigmas = repmat(.25*(tf-ti),N,1);
    
    %% generate fake data
    
    % preallocation
    X = nan(T,N,K);
    R = nan(T,N,K);
    
    % iterate through neurons
    for nn = 1 : N
        for kk = 1 : K
            X(:,nn,kk) = generativerate(...
                t,gammas(nn),mus(nn),lambdas(nn),sigmas(nn));
%             x_padded = generativerate(...
%                 t_padded,gammas(nn),mus(nn),lambdas(nn),sigmas(nn));
%             dur = tf - ti - kernel.paddx(1) + kernel.paddx(2);
%             [n,ts] = poissonprocess(x_padded,dur/t_units);
%             spk_times = ts * t_units + ti + kernel.paddx(1);
%             spk_counts = histcounts(spk_times,'binedges',t_padded);
%             spk_counts = spk_counts / (dt/t_units);
%             R(:,nn,kk) = conv(spk_counts,kernel.pdf,'valid');
            [n,ts] = poissonprocess(X(:,nn,kk),(tf-ti)/t_units);
            spk_times = ts * t_units + ti;
            spk_counts = histcounts(spk_times,'binedges',[t,t(end)+1]);
            R(:,nn,kk) = movsum(spk_counts,dt) / (dt/t_units);
        end
    end
    
    %% compute cross-trial averages
    x = nanmean(X,3);
    r = nanmean(R,3);
    
    %% compute correlation coefficients
    rhos_monotonocity = nan(N,1);
    
    % iterate through neurons
    for nn = 1 : N
        rhos_monotonocity(nn) = corr(t_idcs',r(:,nn));
    end
    
    %% perform linear regression
    mdls = cell(N,1);
    pvals_monotonocity = nan(N,1);
    betas_monotonocity = nan(N,1);
    
    % iterate through neurons
    for nn = 1 : N
        mdls{nn} = fitlm(t_idcs,r(:,nn));
        p = polyfit(t_idcs,r(:,nn)',1);
        pvals_monotonocity(nn) = mdls{nn}.Coefficients.pValue(end);
        betas_monotonocity(nn) = p(1); % mdls{nn}.Coefficients.Estimate(end);
    end
    
    %% stereotypy assessment
    P = 10;
    k = floor(K / P);
    shuffled_trial_idcs = randperm(K,K);
    
    % preallocation
    r_partitions = nan(T,N,P);
    rhos_stereotypy = nan(N,P);
    pvals_stereotypy = nan(N,P);
    
    % iterate through partitions
    for pp = 1 : P
        partition_idcs = shuffled_trial_idcs((1 : k) + (pp - 1) * k);
        r_partitions(:,:,pp) = nanmean(R(:,:,partition_idcs),3);
    end
    
    % compute reference
    r_ref = nanmean(r_partitions,3);
    
    % iterate through partitions
    for pp = 1 : P
        
        % iterate through neurons
        for nn = 1 : N
            [rhos,pvals] = corrcoef(r_ref(:,nn),r_partitions(:,nn,pp));
            rhos_stereotypy(nn,pp) = rhos(1,2);
            pvals_stereotypy(nn,pp) = pvals(1,2);
        end
    end
    
    %     for nn = 1 : N
    %         figure('windowstyle','docked');
    %         hold on;
    %         plot(t,r_ref(:,nn),'linewidth',2);
    %         for pp = 1 : P
    %             plot(t,r_partitions(:,nn,pp),'linewidth',1);
    %         end
    %         text(.05,.95,...
    %             sprintf('\\gamma=%.1f\n\\lambda=%.1f',...
    %             gammas(nn),lambdas(nn)),...
    %             'units','normalized')
    %     end
    
    %% neuron selection
    
    % flag "monotonic" neurons
    monotonicity_flags = ...
        abs(rhos_monotonocity) > rho_monotonocity_cutoff & ...
        abs(betas_monotonocity) > beta_monotonocity_cutoff & ...
        pvals_monotonocity <= pval_monotonocity_cutoff;
    
    % flag stereotypical neurons
    stereotypy_flags = ...
        mean(rhos_stereotypy > rho_stereotypy_cutoff,2) == 1 & ...
        mean(pvals_stereotypy < pval_stereotypy_cutoff,2) == 1;
    
    % flag "ramping" neurons
    ramp_flags = ...
        monotonicity_flags & ...
        stereotypy_flags;
    
    [mean(monotonicity_flags),mean(stereotypy_flags),mean(ramp_flags)]
    
    % store proportion of ramping neurons
    p_ramp(bb) = nanmean(ramp_flags);
    
    %% naive bayes decoder
    %     tensor_ramp = cat(3,r(:,ramp_flags),r(:,ramp_flags));
    %     tensor_non = cat(3,r(:,~ramp_flags),r(:,~ramp_flags));
    train_flags = ismember(1:K,randperm(K,round(K/2)));
    tensor_ramp = cat(3,...
        nanmean(R(:,ramp_flags,train_flags),3),...
        nanmean(R(:,ramp_flags,~train_flags),3));
    tensor_non = cat(3,...
        nanmean(R(:,~ramp_flags,train_flags),3),...
        nanmean(R(:,~ramp_flags,~train_flags),3));
    
    opt = struct();
    opt.n_xpoints = 100;
    opt.time = t;
    opt.train.trial_idcs = 1;
    opt.train.n_trials = numel(opt.train.trial_idcs);
    opt.test.trial_idcs = 2;
    opt.test.n_trials = numel(opt.test.trial_idcs);
    opt.assumepoissonmdl = true;
    opt.verbose = false;
    
    % preallocation
    P_Rt = nan(T,N,opt.n_xpoints);
    [P_tR_ramp,P_Rt(:,ramp_flags,:),~,~] = ...
        naivebayestimedecoder(tensor_ramp,opt);
    [P_tR_non,P_Rt(:,~ramp_flags,:),~,~] = ...
        naivebayestimedecoder(tensor_non,opt);
    
    % normalization
    P_tR_ramp = P_tR_ramp ./ nansum(P_tR_ramp,1);
    P_tR_non = P_tR_non ./ nansum(P_tR_non,1);
    
    %% plot encoding models
    % figure;
    % set(gca,...
    %     'xlim',[ti,tf],...
    %     'xdir','normal',...
    %     'ydir','normal',...
    %     'nextplot','add',...
    %     'plotboxaspectratio',[1,1,1]);
    % xlabel('Time (a.u.)');
    % ylabel('Rate (a.u.)');
    %
    % % iterate through units
    % for nn = 1 : N
    %     cla;
    %     r_bounds = neurons(nn).x_bounds;
    %     r_bw = neurons(nn).x_bw;
    %     if range(r_bounds) == 0
    %         continue;
    %     end
    %     ylim(r_bounds);
    %     title(sprintf('Neuron: %i, bw: %.2f',nn,r_bw));
    %     p_Rt = squeeze(P_Rt(:,nn,:));
    %     p_Rt(isnan(p_Rt)) = max(p_Rt(:));
    %     imagesc([ti,tf],r_bounds,p_Rt');
    %     plot(t,neurons(nn).x_mu,...
    %         'color','w',...
    %         'linewidth',1);
    %     drawnow;
    %     pause(.1);
    % end
    % close;
    
    %% compute decoding error
    
    % preallocation
    mu_ramp = nan(T,1);
    mu_non = nan(T,1);
    sd_ramp = nan(T,1);
    sd_non = nan(T,1);
    
    % iterate through time points
    for tt = 1 : T
        mu_ramp(tt) = P_tR_ramp(:,tt)' * t';
        mu_non(tt) = P_tR_non(:,tt)' * t';
        sd_ramp(tt) = sqrt(P_tR_ramp(:,tt)' * (mu_ramp(tt) - t') .^ 2);
        sd_non(tt) = sqrt(P_tR_non(:,tt)' * (mu_non(tt) - t') .^ 2);
    end
    
    % store current bootstrap
    MU_ramp(:,bb) = mu_ramp;
    MU_non(:,bb) = mu_non;
    SD_ramp(bb) = nanmean(sd_ramp);
    SD_non(bb) = nanmean(sd_non);
    GAMMAS(:,bb) = gammas;
    LAMBDAS(:,bb) = lambdas;
    RAMP_FLAGS(:,bb) = ramp_flags;
    
    %%
    if bb < B
        continue;
    end
    
    %% plot single-neuron averages
    figure('position',[40.2000 41.8000 217.6000 740.8000]);
    set(gca,...
        'xlim',[ti,tf],...
        'ycolor','none',...
        'nextplot','add');
    xlabel('time (a.u.)');
    
    % plot average firing rate
    plot(t,x+(1:N)*quantile(gammas,.85),...
        'color','k',...
        'linewidth',1);
    plot(t,r+(1:N)*quantile(gammas,.85),...
        'color','r',...
        'linewidth',1);
    
    %% tiling
    
    % figure initialization
    figure(figopt,...
        'name','tiling');
    
    % axes initialization
    set(gca,...
        axesopt.default,...
        'xlim',[ti,tf],...
        'ylim',[1,N],...
        'ytick',[1,N],...
        'colormap',hot(2^8),...
        'layer','top',...
        'tickdir','out',...
        'nextplot','add',...
        'plotboxaspectratio',[1,1,1],...
        'linewidth',2,...
        'fontsize',12,...
        'ticklength',[1,1]*.025);
    xlabel('Time (s)');
    ylabel('Neuron #');
    
    % color limits
    clim = [-2,4];
    
    % normalization
    z = zscore(x);
    
    % plot psth raster
    imagesc(t,[1,N],z',clim);
    %     imagesc(t,[1,N],x');
    % imagesc(t,[1,N],x_non(:,:,2)');
    
    % color bar
    clrbar = colorbar;
    clrbar.Ticks = unique([0,clim]);
    clrlabel.string = 'Firing rate (z-score)';
    clrlabel.fontsize = axesopt.default.fontsize * 1.1;
    clrlabel.rotation = 270;
    clrlabel.position = [4,sum(clrbar.Limits)/2,0];
    set(clrbar,...
        axesopt.colorbar,...
        'color','k',...
        'fontsize',axesopt.default.fontsize);
    set(clrbar.Label,...
        'color','k',...
        clrlabel);
    
    %% plot single-neuron averages
    figure(figopt,...
        'position',[300,350,750,420],...
        'color','w');
    set(gca,...
        axesopt.default,...
        'xlim',[ti,tf],...
        'ylim',[0,max(gammas)] + [-1,1] * 0.05 * max(gammas),...
        'ycolor','none',...
        'layer','top',...
        'tickdir','out',...
        'nextplot','add',...
        'plotboxaspectratio',[1,1,1],...
        'linewidth',2,...
        'fontsize',12,...
        'ticklength',[1,1]*.025);
    xlabel('Time (ms)');
    
    % iterate through neurons
    for nn = 1 : N
        
        % plot single-neuron average
        plot(t,x(:,nn),...
            'color',clrs(nn,:),...
            'linewidth',1);
        
        % plot model prediction
        %         plot(t,mdls{nn}.predict(t_idcs'),...
        %             'color',clrs(nn,:),...
        %             'linestyle',repmat('-',1+(pvals_monotonocity(nn)>pval_monotonocity_cutoff),1),...
        %             'linewidth',1);
        
        % annotate correlation coefficient
        text(1.05,nn/N,sprintf('\\rho = %.2f',rhos_monotonocity(nn)),...
            'color',clrs(nn,:),...
            'horizontalalignment','left',...
            'units','normalized');
        
        % annotate regression statistics
        text(-.05,nn/N,...
            sprintf('\\beta = %.2f, p-value = %.2f',betas_monotonocity(nn),pvals_monotonocity(nn)),...
            'color',clrs(nn,:),...
            'horizontalalignment','right',...
            'units','normalized');
        
        % afford ramping class
        if ramp_flags(nn)
            plot(mus(nn),max(ylim),...
                'marker','v',...
                'markerfacecolor',clrs(nn,:),...
                'markeredgecolor','w',....
                'markersize',7.5,...
                'linewidth',1);
        end
    end
    
    %% parameter scatter
    
    % figure initialization
    fig = figure(figopt,...
        'name','scatter');
    
    % axes initialization
    set(gca,axesopt.default,...
        'xscale','log',...
        'yscale','log',...
        'xlimspec','tight',...
        'ylimspec','tight',...
        'zlimspec','tight',...
        'layer','bottom');
    xlabel('log(\gamma)');
    ylabel('log(\lambda)');
    %     xlabel('\gamma');
    %     ylabel('\lambda');
    
    % scatter
    all_gammas = GAMMAS(:);
    all_lambdas = LAMBDAS(:);
    all_flags = RAMP_FLAGS(:) == 1;
    grapeplot(all_gammas,all_lambdas,...
        'markerfacecolor',[1,1,1]);
    grapeplot(all_gammas(all_flags),all_lambdas(all_flags),...
        'markerfacecolor',[1,0,0]);
    
    %     scatterhist(all_gammas,all_lambdas,...
    %         'group',categorical(all_flags,[0,1],{'Non-ramping','Ramping'}),...
    %         'kernel','on',...
    %         'linewidth',1.5);
    
    % legend
    legend({'Non-ramping','Ramping'},...
        'location','southwest',...
        'box','off');
    
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
    
    %% plot posterior averages
    figure(...
        'name','condition-split posterior averages',...
        'numbertitle','off',...
        'windowstyle','docked');
    n_rows = 1;
    n_cols = 2;
    sps = gobjects(n_rows,n_cols);
    for rr = 1 : n_rows
        for cc = 1 : n_cols
            sp_idx = cc + (rr - 1) * n_cols;
            sps(rr,cc) = subplot(n_rows,n_cols,sp_idx);
            xlabel(sps(rr,cc),'real time (a.u.)');
            ylabel(sps(rr,cc),'decoded time (a.u.)');
        end
    end
    set(sps,...
        'xlim',[ti,tf],...
        'ylim',[ti,tf],...
        'xdir','normal',...
        'ydir','normal',...
        'nextplot','add',...
        'plotboxaspectratio',[1,1,1]);
    linkaxes(sps);
    
    % iterate through conditions
    clims = quantile([P_tR_ramp(:);P_tR_non(:)],[0,1]);
    
    % ramping posteriors
    title(sps(1),sprintf('Ramping neurons (%i/%i)',sum(ramp_flags),N));
    imagesc(sps(1),[ti,tf],[ti,tf],P_tR_ramp',clims);
    plot(sps(1),xlim(sps(1)),ylim(sps(1)),'-k');
    plot(sps(1),xlim(sps(1)),ylim(sps(1)),'--w');
    
    % non-ramping posteriors
    title(sps(2),sprintf('Non-ramping neurons (%i/%i)',sum(~ramp_flags),N));
    imagesc(sps(2),[ti,tf],[ti,tf],P_tR_non',clims);
    plot(sps(2),xlim(sps(2)),ylim(sps(2)),'-k');
    plot(sps(2),xlim(sps(2)),ylim(sps(2)),'--w');
    
    %%
    figure(...
        'name','posterior_means_SD',...
        'numbertitle','off',...
        'windowstyle','docked');
    n_rows = 1;
    n_cols = 2;
    sps = gobjects(n_rows,n_cols);
    for rr = 1 : n_rows
        for cc = 1 : n_cols
            sp_idx = cc + (rr - 1) * n_cols;
            sps(rr,cc) = subplot(n_rows,n_cols,sp_idx);
            xlabel(sps(rr,cc),'Real time (a.u.)');
            ylabel(sps(rr,cc),'Decoded time (a.u.)');
        end
    end
    set(sps,...
        'xlim',[ti,tf],...
        'ylim',[ti,tf],...
        'xdir','normal',...
        'ydir','normal',...
        'nextplot','add',...
        'plotboxaspectratio',[1,1,1]);
    linkaxes(sps);
    
    % ramping posteriors
    title(sps(1),sprintf('Ramping neurons (%i/%i)',sum(ramp_flags),N));
    errorbar(sps(1),t,mu_ramp,sd_ramp,...
        'color','r',...
        'capsize',0);
    % plot(sps(1),xlim(sps(1)),ylim(sps(1)),'--k');
    
    % non-ramping posteriors
    title(sps(2),sprintf('Non-ramping neurons (%i/%i)',sum(~ramp_flags),N));
    errorbar(sps(2),t,mu_non,sd_non,...
        'color','k',...
        'capsize',0);
    % plot(sps(2),xlim(sps(2)),ylim(sps(2)),'--k');

    %%
    
    figure;
    hold on;
    xlim([ti,tf]);
    
    for tt = 1 : T
        cla;
        plot(t,P_tR_ramp(:,tt),...
            'color','r',...
            'linewidth',1.5);
        plot(mu_ramp(tt),max(P_tR_ramp(:,tt)),...
            'color','r',...
            'marker','v',...
            'markersize',7.5,...
            'linewidth',1.5);
        plot(mu_ramp(tt)+[-1,1]*sd_ramp(tt),[1,1]*max(P_tR_ramp(:,tt)),...
            'color','r',...
            'linewidth',1.5);
        
        plot(t,P_tR_non(:,tt),...
            'color','k',...
            'linewidth',1.5);
        plot(mu_non(tt),max(P_tR_non(:,tt)),...
            'color','k',...
            'marker','v',...
            'markersize',7.5,...
            'linewidth',1.5);
        plot(mu_non(tt)+[-1,1]*sd_non(tt),[1,1]*max(P_tR_non(:,tt)),...
            'color','k',...
            'linewidth',1.5);
        
        if B == B
            pause(.05)
        end
    end
    %%
    % figure;
    % hold on;
    %
    % x_test = linspace(0,1,T);
    % x_pdf = normpdf(x_test,.35,.05);
    % x_pdf = x_pdf / nansum(x_pdf);
    %
    % x_mu = x_pdf * x_test';
    % x_sd = sqrt(x_pdf * (x_mu - x_test') .^ 2);
    %
    % plot(x_test,x_pdf);
    % plot(x_mu,max(x_pdf),...
    %     'color','r',...
    %     'marker','v',...
    %     'markersize',7.5,...
    %     'linewidth',1.5)
    % plot(x_mu+[-1,1]*x_sd,[1,1]*max(x_pdf),...
    %     'color','r',...
    %     'linewidth',1.5)
    
    %%
    
    figure;
    hold on;
    xlabel('SD P(t | R_{ramping})');
    ylabel('SD P(t | R_{non-ramping})');
    
    plot(sd_ramp,sd_non,'.');
    plot(sd_ramp(1),sd_non(1),'ok');
    plot(nanmean(sd_ramp),nanmean(sd_non),'xr');
    
    % [~,idx] = max(P_tR_non,[],2);
    % plot(t,t(idx),'linewidth',2);
    
    axis tight;
    axis square;
    lims = [min([xlim,ylim],[],'all'),max([xlim,ylim],[],'all')];
    plot(lims,lims,'--k')
end

%%
figure;
hold on;
xlabel('SD (ms)');
[~,edges] = histcounts([SD_ramp,SD_non]);
histogram(SD_ramp,edges,...
    'facecolor','k');
histogram(SD_non,edges,...
    'facecolor','r');

%%
figure(figopt);
axes(axesopt.default,...
    'xlim',[-.5,1.5],...
    'xtick',[0,1],...
    'xticklabel',{'Ramping','Non-ramping'},...
    'plotboxaspectratio',[1,2,1]);
ylabel('SD (ms)');
errorbar(0,nanmean(SD_ramp),nanstd(SD_ramp),...
    'color','r',...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor','r',...
    'linewidth',1.5);
errorbar(1,nanmean(SD_non),nanstd(SD_non),...
    'color','k',...
    'marker','o',...
    'markersize',7.5,...
    'markeredgecolor','k',...
    'markerfacecolor','k',...
    'linewidth',1.5);

%%
figure(figopt,...
    'name','posterior_means',...
    'numbertitle','off',...
    'windowstyle','docked');
n_rows = 1;
n_cols = 2;
sps = gobjects(n_rows,n_cols);
for rr = 1 : n_rows
    for cc = 1 : n_cols
        sp_idx = cc + (rr - 1) * n_cols;
        sps(rr,cc) = subplot(n_rows,n_cols,sp_idx);
        xlabel(sps(rr,cc),'Real time (ms)');
        ylabel(sps(rr,cc),'Decoded time (ms)');
    end
end
set(sps,...
    axesopt.default,...
    'xlim',[ti,tf],...
    'ylim',[ti,tf],...
    'xdir','normal',...
    'ydir','normal',...
    'nextplot','add',...
    'plotboxaspectratio',[1,1,1]);
linkaxes(sps);

% ramping posteriors
title(sps(1),sprintf('Ramping neurons (%.0f%%)',nanmean(p_ramp)*100));
errorbar(sps(1),t,nanmean(MU_ramp,2),nanstd(MU_ramp,0,2),...
    'color',rgb('crimson'),...
    'linewidth',.1,...
    'capsize',0);
plot(sps(1),t,nanmean(MU_ramp,2),...
    'color','k',...
    'linewidth',2);
% plot(sps(1),xlim(sps(1)),ylim(sps(1)),'--k');

% non-ramping posteriors
title(sps(2),sprintf('Non-ramping neurons (%.0f%%)',(1-nanmean(p_ramp))*100));
errorbar(sps(2),t,nanmean(MU_non,2),nanstd(MU_non,0,2),...
    'color',rgb('crimson'),...
    'linewidth',.1,...
    'capsize',0);
plot(sps(2),t,nanmean(MU_non,2),...
    'color','k',...
    'linewidth',2);
% plot(sps(2),xlim(sps(2)),ylim(sps(2)),'--k');

%%
figure;
axes(axesopt.default);
xlabel('Firing rate (Hz)');
histogram(nanmean(R,[1,3]),round(N/2),...
    'facealpha',1,...
    'edgecolor','w',...
    'facecolor','k',...
    'linewidth',1.5);

%% firing rate function
function x = generativerate(time,gamma,mu,lambda,sigma)
% mu_pd = truncate(...
%     makedist('normal','mu',mu,'sigma',lambda),time(1),time(end));
% if lambda == 0
%     mu_sample = mu;
% else
%     mu_sample = random(mu_pd);
% end
mu_sample = clamp(normrnd(mu,lambda),time(1),time(end));
x_pdf = normpdf(time,mu_sample,sigma);
x = gamma * x_pdf ./ max(x_pdf);
end