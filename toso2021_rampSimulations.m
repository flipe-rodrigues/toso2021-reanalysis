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
% T = 100; 	% time points
N = 50;     % neurons
K = 100;    % trials
B = 100;   	% bootstrap iterations
P = 4;      % number of partitions with which to assess stereotypy

%% time settings
dt = 100;
ti = 0;
tf = 1000;
t = ti : tf;
t_padded = ti - dt/2  : tf + dt/2;
T = numel(t);
% t = linspace(ti,tf,T);
t_idcs = 1 : T;
t_units = 1e3;
% dt = diff(t(1:2));

%% "ramping" criteria

% "monotonocity" criteria
rho_monotonocity_cutoff = .5;
% beta_monotonocity_cutoff = .05;     % what they say they do in the methods
beta_monotonocity_cutoff = .004;	% what's in their code
pval_monotonocity_cutoff = .05;

% stereotypy criteria
rho_stereotypy_cutoff = .5;
pval_stereotypy_cutoff = .01;

%% model settings

% model parameters
gamma_ranges = [...
    [5,5]; ...
    [5,5]; ...
    [0,100]; ...
    [0,100]];
lambda_ranges = [...
    [0,0]; ...
    [0,tf]; ...
    [0,0]; ...
    [0,tf]];
M = size(gamma_ranges,1);

% model labels
model_labels = cell(M,1);
for mm = 1 : M
    model_labels{mm} = sprintf([...
        '\\gamma\\in[%i,%i]','\\newline',...
        '\\lambda\\in[%i,%i]'],...
        gamma_ranges(mm,1),gamma_ranges(mm,2),...
        lambda_ranges(mm,1),lambda_ranges(mm,2));
end
disp(model_labels);

%% model simulations

% preallocation
MU_ramp = nan(T,B,M);
MU_non = nan(T,B,M);
SD_ramp = nan(T,B,M);
SD_non = nan(T,B,M);
P_RAMP = nan(B,M);
RAMP_FLAGS = nan(N,B,M);
MUS = nan(N,B,M);
GAMMAS = nan(N,B,M);
LAMBDAS = nan(N,B,M);

close all

% iterate through models
for mm = 1 : M
    
    % iterate through bootstrap iterations
    for bb = 1 : B
        progressreport(bb,B,'bootstrapping')
        
        %% model parameters
        gammas = clamp(exprnd(10,N,1),gamma_ranges(mm,1),gamma_ranges(mm,2));
        lambdas = clamp(exprnd(.15*(tf-ti),N,1),lambda_ranges(mm,1),lambda_ranges(mm,2));
        mus = linspace(ti,tf,N)';
%         sigmas = repmat(.25*(tf-ti),N,1);
        sigmas = abs(normrnd(.25,.05,N,1)) * (tf - ti);
        
        %% generate fake data
        
        % preallocation
        X = nan(T,N,K);
        R = nan(T,N,K);
        
        % baseline firing rate
        bsl_fr = 2;
        
        % iterate through neurons
        for nn = 1 : N
            for kk = 1 : K
                X(:,nn,kk) = bsl_fr + generativerate(...
                    t,gammas(nn),mus(nn),lambdas(nn),sigmas(nn));
                x_padded = bsl_fr + generativerate(...
                    t_padded,gammas(nn),mus(nn),lambdas(nn),sigmas(nn));
                dur_padded = range(t_padded) / t_units;
                [~,ts] = poissonprocess(x_padded,dur_padded);
                spk_times = ts * t_units + t_padded(1);
                spk_counts = histcounts(spk_times,'binedges',t_padded);
                r_padded = movsum(spk_counts,dt) / (dt / t_units);
                t_flags = t_padded >= ti & t_padded <= tf;
                R(:,nn,kk) = r_padded(t_flags);
            end
        end
        
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        R = X;
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
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
            betas_monotonocity(nn) = p(1);
        end
        
        %% stereotypy assessment
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
        
        %% naive bayes decoder
        train_flags = ismember(1:K,randperm(K,round(K/2)));
        tensor_ramp = cat(3,...
            nanmean(R(:,ramp_flags,train_flags),3),...
            nanmean(R(:,ramp_flags,~train_flags),3));
        tensor_non = cat(3,...
            nanmean(R(:,~ramp_flags,train_flags),3),...
            nanmean(R(:,~ramp_flags,~train_flags),3));
        
        % enforce equal numbers of neurons on both clusters
        n_ramp = sum(ramp_flags);
        n_non = sum(~ramp_flags);
        n = min(n_ramp,n_non);
        tensor_ramp = tensor_ramp(:,randperm(n_ramp,n),:);
        tensor_non = tensor_non(:,randperm(n_non,n),:);

        % decoding options
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
        P_tR_ramp = naivebayestimedecoder(tensor_ramp,opt);
        P_tR_non = naivebayestimedecoder(tensor_non,opt);
        
        % normalization
        P_tR_ramp = P_tR_ramp ./ nansum(P_tR_ramp,1);
        P_tR_non = P_tR_non ./ nansum(P_tR_non,1);
        
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
        
        %% store current bootstrap run
        MU_ramp(:,bb,mm) = mu_ramp;
        MU_non(:,bb,mm) = mu_non;
        SD_ramp(:,bb,mm) = sd_ramp;
        SD_non(:,bb,mm) = sd_non;
        P_RAMP(bb,mm) = nanmean(ramp_flags);
        MUS(:,bb,mm) = mus;
        GAMMAS(:,bb,mm) = gammas;
        LAMBDAS(:,bb,mm) = lambdas;
        RAMP_FLAGS(:,bb,mm) = ramp_flags;
        
        %% skip plotting single-run stuff if it's not the last run
        if bb < B
            continue;
        end
        
        %% plot single-neuron averages
        figure(figopt,...
            'position',[40.2000 41.8000 217.6000 740.8000],...
            'name',sprintf('generative_vs_measured_%i',mm));
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
            'position',[100,343,1325,420],...
            'name',sprintf('tiling_%i',mm));
        
        % axes initialization
        n_sps = 3;
        sps = gobjects(n_sps,1);
        for ii = 1 : n_sps
            sps(ii) = subplot(1,n_sps,ii);
            xlabel(sps(ii),'Time (s)');
            ylabel(sps(ii),'Neuron #');
        end
        set(sps,...
            axesopt.default,...
            'xlim',[ti,tf],...
            'colormap',hot(2^8),...
            'layer','top',...
            'tickdir','out',...
            'nextplot','add',...
            'plotboxaspectratio',[1,1,1],...
            'linewidth',2,...
            'fontsize',12,...
            'ticklength',[1,1]*.025);
        set(sps(1),...
            'ylim',[1,N],...
            'ytick',[1,N]);
        set(sps(2),...
            'ylim',[1,sum(ramp_flags)],...
            'ytick',[1,sum(ramp_flags)]);
        set(sps(3),...
            'ylim',[1,sum(~ramp_flags)],...
            'ytick',[1,sum(~ramp_flags)]);
        title(sps(1),'All');
        title(sps(2),'Ramps');
        title(sps(3),'Non-ramps');
        
        % normalization
        z = zscore(x);
        
        % color limits
        r_clim = [min(r,[],'all'),max(r,[],'all')];
        z_clim = [-1,1] .* max(abs(z),[],'all');
        
        % plot psth raster
        imagesc(sps(1),t,[1,N],r',r_clim);
        imagesc(sps(2),t,[1,sum(ramp_flags)],r(:,ramp_flags)',r_clim);
        imagesc(sps(3),t,[1,sum(~ramp_flags)],r(:,~ramp_flags)',r_clim);
%         imagesc(sps(2),t,[1,N],z(:,ramp_flags)',clim);
%         imagesc(sps(3),t,[1,N],z(:,~ramp_flags)',clim);
        
        % color bar
%         clrbar = colorbar;
% %         clrbar.Ticks = unique([0,clim]);
% %         clrlabel.string = 'Firing rate (z-score)';
%         clrlabel.string = 'Firing rate (Hz)';
%         clrlabel.fontsize = axesopt.default.fontsize * 1.1;
%         clrlabel.rotation = 270;
%         clrlabel.position = [4,sum(clrbar.Limits)/2,0];
%         set(clrbar,...
%             axesopt.colorbar,...
%             'color','k',...
%             'fontsize',axesopt.default.fontsize);
%         set(clrbar.Label,...
%             'color','k',...
%             clrlabel);
        
        %% plot single-neuron averages
        figure(figopt,...
            'position',[300,350,750,420],...
            'name',sprintf('single_neurons_%i',mm));
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
        xlabel('Time (a.u.)');
        
        % color settings
        clrs = cool(N);
        
        % iterate through neurons
        for nn = 1 : N
            
            % plot single-neuron average
            plot(t,x(:,nn),...
                'color',clrs(nn,:),...
                'linewidth',1);
            
            % plot regression prediction
            %             plot(t,mdls{nn}.predict(t_idcs'),...
            %                 'color',clrs(nn,:),...
            %                 'linestyle',repmat('-',1+(pvals_monotonocity(nn)>pval_monotonocity_cutoff),1),...
            %                 'linewidth',1);
            
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
        figure(figopt,...
            'name',sprintf('parameter_scatter_%i',mm));
        
        % axes initialization
        set(gca,axesopt.default,...
            'xscale','linear',...
            'yscale','linear',...
            ...'xscale','log',...
            ...'yscale','log',...
            'xlimspec','tight',...
            'ylimspec','tight',...
            'zlimspec','tight',...
            'layer','bottom');
        %         xlabel('log(\gamma)');
        %         ylabel('log(\lambda)');
        xlabel('\gamma');
        ylabel('\lambda');
        
        % scatter
        all_gammas = GAMMAS(:,:,mm);
        all_gammas = all_gammas(:);
        all_lambdas = LAMBDAS(:,:,mm);
        all_lambdas = all_lambdas(:);
        all_flags = RAMP_FLAGS(:,:,mm) == 1;
        all_flags = all_flags(:);
        grapeplot(all_gammas,all_lambdas,...
            'markerfacecolor',[1,1,1]);
        grapeplot(all_gammas(all_flags),all_lambdas(all_flags),...
            'markerfacecolor',[1,0,0]);
        
        % legend
        legend({'Non-ramping','Ramping'},...
            'location','best',...
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
            'name',sprintf('posterior_averages_%i',mm),...
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
        clim = [0,1/T*N]; % quantile([P_tR_ramp(:);P_tR_non(:)],[0,1]);
        
        % ramping posteriors
        title(sps(1),sprintf('Ramping neurons (%i/%i)',sum(ramp_flags),N));
        imagesc(sps(1),[ti,tf],[ti,tf],P_tR_ramp',clim);
        plot(sps(1),xlim(sps(1)),ylim(sps(1)),'-k');
        plot(sps(1),xlim(sps(1)),ylim(sps(1)),'--w');
        
        % non-ramping posteriors
        title(sps(2),sprintf('Non-ramping neurons (%i/%i)',sum(~ramp_flags),N));
        imagesc(sps(2),[ti,tf],[ti,tf],P_tR_non',clim);
        plot(sps(2),xlim(sps(2)),ylim(sps(2)),'-k');
        plot(sps(2),xlim(sps(2)),ylim(sps(2)),'--w');
        
        %%
        figure(...
            'name',sprintf('posterior_means_SD_%i',mm),...
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
            'color',ramp_clrs(1,:),...
            'capsize',0);
        % plot(sps(1),xlim(sps(1)),ylim(sps(1)),'--k');
        
        % non-ramping posteriors
        title(sps(2),sprintf('Non-ramping neurons (%i/%i)',sum(~ramp_flags),N));
        errorbar(sps(2),t,mu_non,sd_non,...
            'color',ramp_clrs(2,:),...
            'capsize',0);
        % plot(sps(2),xlim(sps(2)),ylim(sps(2)),'--k');
        
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
end

%% decoding error across models

% figure initialization
fig = figure(figopt,...
    'name','ramps_sd_models');

% axes initialization
xxtick = unique((1:M)+[-1;0;1]*.05*M);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:M)) = {''};
xxticklabel(ismember(xxtick,1:M)) = model_labels;
axes(axesopt.default,...
    'plotboxaspectratio',[1,1,1],...
    'color','none',...
    'xlim',[1,M]+[-1,1]*.1*M,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'clipping','off',...
    'layer','bottom');
xlabel('Parameter range');
ylabel('Decoding precision SD (a.u.)');

% offset between ramps and non-ramps
xoffsets = [-1,1] * .15;

% choice of precision function
precisionfun = @(x,d) 1 ./ nanvar(x,0,d);
precisionfun = @(x,d) nanstd(x,0,d);

% choice of average and error functions
avgfun = @(x) nanmean(x);
errfun = @(x) [1,1] .* nanstd(x);
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);

%
X_ramp = MU_ramp;
X_non = MU_non;
X_ramp = SD_ramp;
X_non = SD_non;

% iterate through models
for mm = 1 : M
    sd_ramp_boots = precisionfun(X_ramp(:,:,mm),2);
    sd_non_boots = precisionfun(X_non(:,:,mm),2);
    sd_ramp_avg = avgfun(sd_ramp_boots);
    sd_non_avg = avgfun(sd_non_boots);
    sd_ramp_err = errfun(sd_ramp_boots);
    sd_non_err = errfun(sd_non_boots);
    plot(mm+xoffsets,[sd_ramp_avg,sd_non_avg],...
        'color','k',...
        'linewidth',1.5);
    errorbar(mm+xoffsets(1),sd_ramp_avg,...
        sd_ramp_err(1),sd_ramp_err(2),...
        'color','k',...
        'marker','o',...
        'markersize',7.5,...
        'markeredgecolor','k',...
        'markerfacecolor',ramp_clrs(1,:),...
        'linewidth',1.5,...
        'capsize',0);
    errorbar(mm+xoffsets(2),sd_non_avg,...
        sd_non_err(1),sd_non_err(2),...
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
yytick = linspace(yylim(1),yylim(2),5);
yyticklabel = num2cell(yytick);
set(gca,...
    'ylim',yylim + [-1,1] * .05 * range(yylim),...
    'ytick',yytick,...
    'yticklabel',yyticklabel);

% iterate through alignments
for mm = 1 : M
    sd_ramp_boots = precisionfun(X_ramp(:,:,mm),2);
    sd_non_boots = precisionfun(X_non(:,:,mm),2);
    xx = [-1,1] * .5 / 3 + mm;
    yy = [1,1] * max(yylim);
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    [~,pval] = kstest2(sd_ramp_boots,sd_non_boots);
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
            xlabel(sps(rr,cc),'Time (a.u.)');
            ylabel(sps(rr,cc),'Decoded time (a.u.)');
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
    title(sps(1),sprintf('Ramping neurons (%.0f%%)',nanmean(P_RAMP(:,mm))*100));
    errorbar(sps(1),t,nanmean(MU_ramp(:,:,mm),2),nanstd(MU_ramp(:,:,mm),0,2),...
        'color',ramp_clrs(1,:),...
        'linewidth',.1,...
        'capsize',0);
    plot(sps(1),t,nanmean(MU_ramp(:,:,mm),2),...
        'color','w',...
        'linewidth',1.5);
    for bb = 1 : B
        plot(sps(1),t,nanmean(MU_ramp(:,bb,mm),2),...
            'color','b',...
            'linewidth',.1);
    end
    
    % non-ramping posteriors
    title(sps(2),sprintf('Non-ramping neurons (%.0f%%)',(1-nanmean(P_RAMP(:,mm)))*100));
    errorbar(sps(2),t,nanmean(MU_non(:,:,mm),2),nanstd(MU_non(:,:,mm),0,2),...
        'color',ramp_clrs(2,:),...
        'linewidth',.1,...
        'capsize',0);
    plot(sps(2),t,nanmean(MU_non(:,:,mm),2),...
        'color','w',...
        'linewidth',1.5);
    for bb = 1 : B
        plot(sps(2),t,nanmean(MU_non(:,bb,mm),2),...
            'color','b',...
            'linewidth',.1);
    end
    
    text(sps(1),.05,.95,model_labels{mm},...
        'horizontalalignment','left',...
        'verticalalignment','top',...
        'units','normalized');
    
    % identity lines
    plot(sps(1),xlim(sps(1)),ylim(sps(1)),'--b','linewidth',1.5);
    plot(sps(2),xlim(sps(2)),ylim(sps(2)),'--b','linewidth',1.5);
end

%%
figure;
for mm = 1 : M
    subplot(M,1,mm);
    histogram(P_RAMP(:,mm),(0:.025:1.05)-.0125);
    xlim([0,1]);
end
return;

%% example neurons with opposite gammas & lambdas
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
        [n,ts] = poissonprocess(X_eg(:,nn,kk),(tf-ti)/t_units);
        spk_times = ts * t_units + ti;
        spk_counts = histcounts(spk_times,'binedges',[t,t(end)+1]);
        R_eg(:,nn,kk) = movsum(spk_counts,dt) / (dt/t_units);
    end
end

% compute cross-trial mean
x_eg = nanmean(X_eg,3);

% perform linear regression
pvals_monotonocity_eg = nan(2,1);
betas_monotonocity_eg = nan(2,1);

% iterate through neurons
for nn = 1 : 2
    mdl = fitlm(t_idcs,x_eg(:,nn));
    p = polyfit(t_idcs,x_eg(:,nn)',1);
    pvals_monotonocity_eg(nn) = mdl.Coefficients.pValue(end);
    betas_monotonocity_eg(nn) = p(1);
end

% figure initialization
fig = figure(figopt,...
    'position',[744,630,460,420],...
    'name','ramp_firingrate');

% axes initialization
[~,idx] = max(x_eg(:,1));
axes(axesopt.default,...
    'xcolor','k',...
    'xtick',unique([ti,t(idx),tf]),...
    'xticklabel',{num2str(ti),'\mu_0',num2str(tf)},...
    'ylim',[0,5]+[-1,1]*.15*range([0,5]),...
    'ytick',[0,5],...
    'clipping','off',...
    'plotboxaspectratio',[3,1,1]);
xlabel('Time (a.u.)');
ylabel('Firing rate (Hz)');

plot(t,x_eg(:,1),...
    'color','k',...
    'linewidth',1.5);
plot(t,x_eg(:,2),...
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
xlabel('Time (a.u.)');
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
mu_sample = clamp(normrnd(mu,lambda),time(1),time(end));
x_pdf = normpdf(time,mu_sample,sigma);
x = gamma * x_pdf ./ max(x_pdf);
end