%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% close all
close all

%% notes
% beta depends on time units (s or ms), and firing rate;
% rho and beta are redundant;
% stereoptypy is a circular criterion for the point of decodability

%% seed fixing
% rng(0);

%% tensor settings
N = 1e3;        % neurons
N_clus = 100;	% neurons per cluster
K = 300;        % trials
S = 100;        % simulations
P = 4;          % number of partitions with which to assess stereotypy

%% temporal smoothing kernel
gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);

%% ROI settings
roi = [0,t_set(end)];
roi_decoding = [0,t_set(end-2)];
roi_onset = [-500,500] + roi(1);
roi_offset = [-500,500] + roi(2);
T_roi = range(roi);
T_roi_decoding = range(roi_decoding);

%% time settings
ti = roi_onset(1);
tf = roi_offset(2);
ti_padded = ti + gauss_kernel.paddx(1);
tf_padded = tf + gauss_kernel.paddx(2);
T = (tf - ti) / psthbin;
T_padded = (tf_padded - ti_padded) / psthbin;
t = linspace(ti,tf,T);
t_roi = linspace(roi(1),roi(2),T_roi);
t_roi_decoding = linspace(roi_decoding(1),roi_decoding(2),T_roi_decoding);
t_padded = linspace(ti_padded,tf_padded,T_padded);
t_units = 1e3;
dt = diff(t(1:2));
roi_flags = ...
    t >= roi(1) & ...
    t <= roi(2);
roi_decoding_flags = ...
    t >= roi_decoding(1) & ...
    t <= roi_decoding(2);
roi_onset_flags = ...
    t >= roi_onset(1) & ...
    t <= roi_onset(2);
roi_offset_flags = ...
    t >= roi_offset(1) & ...
    t <= roi_offset(2);

%% "ramping" criteria

% "monotonocity" criteria
rho_monotonocity_cutoff = .5;
% beta_monotonocity_cutoff = .05;     % what they say they do in the methods
% beta_monotonocity_cutoff = .004;	% what's in their code
beta_monotonocity_cutoff = .003;
pval_monotonocity_cutoff = .05;

% stereotypy criteria
rho_stereotypy_cutoff = .5;
pval_stereotypy_cutoff = .01;

%% model settings

% model parameters
gamma_ranges = [...
    [1,1,1]*3; ...
    [1,1,1]*3; ...
    [0,100,3]; ...
    [0,100,3]];
eta_ranges = [...
    [0,0,0]; ...
    [0,range(roi)*2,range(roi)*.15]; ...
    [0,0,0]; ...
    [0,range(roi)*2,range(roi)*.15]];
M = size(eta_ranges,1);

% model labels
model_labels = cell(M,1);
for mm = 1 : M
    model_labels{mm} = sprintf([...
        '\\gamma\\in[%i,%i]','\\newline',...
        '\\eta\\in[%i,%i]'],...
        gamma_ranges(mm,1),gamma_ranges(mm,2),...
        eta_ranges(mm,1),eta_ranges(mm,2));
end
model_labels = {...
    '\eta=\eta_0\newline\gamma=\gamma_0';...
    '\eta\sim\newline\gamma=\gamma_0';...
    '\eta=\eta_0\newline\gamma\sim';...
    '\eta\sim\newline\gamma\sim';...
    };
disp(model_labels);

%% model simulations

% preallocation
P_TR_RAMP = nan(T_roi_decoding,T_roi_decoding,S,M);
P_TR_NON = nan(T_roi_decoding,T_roi_decoding,S,M);
P_TR_ALL = nan(T_roi_decoding,T_roi_decoding,S,M);
MAP_ramp = nan(T_roi_decoding,S,M);
MAP_non = nan(T_roi_decoding,S,M);
MAP_all = nan(T_roi_decoding,S,M);
MU_ramp = nan(T_roi_decoding,S,M);
MU_non = nan(T_roi_decoding,S,M);
MU_all = nan(T_roi_decoding,S,M);
MED_ramp = nan(T_roi_decoding,S,M);
MED_non = nan(T_roi_decoding,S,M);
MED_all = nan(T_roi_decoding,S,M);
SD_ramp = nan(T_roi_decoding,S,M);
SD_non = nan(T_roi_decoding,S,M);
SD_all = nan(T_roi_decoding,S,M);
IQR_ramp = nan(T_roi_decoding,S,M);
IQR_non = nan(T_roi_decoding,S,M);
IQR_all = nan(T_roi_decoding,S,M);
MUS = nan(N_clus*2,S,M);
SIGMAS = nan(N_clus*2,S,M);
GAMMAS = nan(N_clus*2,S,M);
ETAS = nan(N_clus*2,S,M);
MEAN_FR = nan(N_clus*2,S,M);
SELECTED_RAMP_FLAGS = false(N_clus*2,S,M);
ALL_RAMP_FLAGS = false(N,S,M);
P_RAMP = nan(S,M);

% analysis metrics
TUNING = nan(N_clus*2,S,M);
FR_RANGE = nan(N_clus*2,S,M);
STEREOTYPY = nan(N_clus*2,S,M);

% all unit indices
all_idcs = 1 : N;

% iterate through models
for mm = 1 : M
    
    % iterate through simulations
    for ss = 1 : S
        progressreport(ss,S,'simulating')
        
        %% model parameters
        mus = sort(unifrnd(ti,tf,N,1));
        gammas = clamp(...
            exprnd(gamma_ranges(mm,3),N,1),...
            gamma_ranges(mm,1),gamma_ranges(mm,2));
        etas = clamp(...
            exprnd(eta_ranges(mm,3),N,1),...
            eta_ranges(mm,1),eta_ranges(mm,2));
        sigmas = ones(N,1) * .25 * range(roi);
        
        %% generate fake data
        
        % preallocation
        X = nan(T,N,K);
        R = nan(T,N,K);
        
        % sample baseline firing rates
        bsl_frs = datasample(fr_min.s2(flagged_neurons),N);
            
        % iterate through neurons
        for nn = 1 : N
%             neuron_flags = data.NeuronNumb == flagged_neurons(nn);
%             trial_flags = ...
%                 valid_flags & ...
%                 neuron_flags;
%             n_trials = sum(trial_flags);
            
            % sample S2 durations
            t2_samples = datasample(t2(valid_flags),K);

            % iterate through trials
            for kk = 1 : K
                X(:,nn,kk) = bsl_frs(nn) + generativerate(...
                    t,gammas(nn),mus(nn),etas(nn),sigmas(nn));
                x_padded = bsl_frs(nn) + generativerate(...
                    t_padded,gammas(nn),mus(nn),etas(nn),sigmas(nn));
                dur_padded = (tf_padded - ti_padded) / t_units;
                [~,ts] = poissonprocess(x_padded,dur_padded);
                spk_times = ts * t_units + ti_padded;
                spk_counts = histcounts(spk_times,'binedges',t_padded);
                R(:,nn,kk) = conv(spk_counts/(dt/t_units),gauss_kernel.pdf,'valid');
                
                % simulate variable stimulus durations
                time_flags = ...
                    t <= t2_samples(kk) | ...
                    t >= roi(end);
                R(~time_flags,nn,kk) = nan;
            end
        end
        
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%         R = X;
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        %% compute cross-trial averages
        x = nanmean(X,3);
        r = nanmean(R,3);
        
        %% compute correlation coefficients
        rhos_monotonocity_onset = nan(N,1);
        rhos_monotonocity_offset = nan(N,1);
        
        % iterate through neurons
        for nn = 1 : N
            rhos_monotonocity_onset(nn) = ...
                corr((1:T_roi)',r(roi_onset_flags,nn));
            rhos_monotonocity_offset(nn) = ...
                corr((1:T_roi)',r(roi_offset_flags,nn));
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
            mdls_onset{nn} = fitlm(1:T_roi,r(roi_onset_flags,nn));
            mdls_offset{nn} = fitlm(1:T_roi,r(roi_offset_flags,nn));
            p_onset = polyfit(1:T_roi,r(roi_onset_flags,nn)',1);
            p_offset = polyfit(1:T_roi,r(roi_offset_flags,nn)',1);
            pvals_monotonocity_onset(nn) = ...
                mdls_onset{nn}.Coefficients.pValue(end);
            pvals_monotonocity_offset(nn) = ...
                mdls_offset{nn}.Coefficients.pValue(end);
            betas_monotonocity_onset(nn) = p_onset(1);
            betas_monotonocity_offset(nn) = p_offset(1);
        end
        
        %% stereotypy assessment
        k = floor(K / P);
        shuffled_trial_idcs = randperm(K,K);
        
        % preallocation
        r_partitions_onset = nan(T_roi,N,P);
        r_partitions_offset = nan(T_roi,N,P);
        rhos_stereotypy_onset = nan(N,P);
        rhos_stereotypy_offset = nan(N,P);
        pvals_stereotypy_onset = nan(N,P);
        pvals_stereotypy_offset = nan(N,P);
        
        % iterate through partitions
        for pp = 1 : P
            partition_idcs = shuffled_trial_idcs((1 : k) + (pp - 1) * k);
            r_partitions_onset(:,:,pp) = ...
                nanmean(R(roi_onset_flags,:,partition_idcs),3);
            r_partitions_offset(:,:,pp) = ...
                nanmean(R(roi_offset_flags,:,partition_idcs),3);
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
        train_flags = ismember(1:K,randperm(K,round(K/2)));
        tensor_ramp = cat(3,...
            nanmean(R(roi_decoding_flags,ramp_flags,train_flags),3),...
            nanmean(R(roi_decoding_flags,ramp_flags,~train_flags),3));
        tensor_non = cat(3,...
            nanmean(R(roi_decoding_flags,~ramp_flags,train_flags),3),...
            nanmean(R(roi_decoding_flags,~ramp_flags,~train_flags),3));
%         tensor_all = cat(3,...
%             nanmean(R(roi_flags,:,train_flags),3),...
%             nanmean(R(roi_flags,:,~train_flags),3));
        
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
%         tensor_all = tensor_all(:,selected_idcs(1:N_clus),:);
        
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
%         P_tR_all = naivebayestimedecoder(tensor_all,opt);
        
        %% compute decoding statistics
        
        % preallocation
        mu_ramp = nan(T_roi_decoding,1);
        mu_non = nan(T_roi_decoding,1);
%         mu_all = nan(T_roi_decoding,1);
        sd_ramp = nan(T_roi_decoding,1);
        sd_non = nan(T_roi_decoding,1);
%         sd_all = nan(T_roi_decoding,1);
        
        % iterate through time points
        for tt = 1 : T_roi_decoding
            mu_ramp(tt) = P_tR_ramp(tt,:) * t_roi_decoding';
            mu_non(tt) = P_tR_non(tt,:) * t_roi_decoding';
%             mu_all(tt) = P_tR_all(tt,:) * t_roi_decoding';
            sd_ramp(tt) = sqrt(P_tR_ramp(tt,:) * (mu_ramp(tt) - t_roi_decoding') .^ 2);
            sd_non(tt) = sqrt(P_tR_non(tt,:) * (mu_non(tt) - t_roi_decoding') .^ 2);
%             sd_all(tt) = sqrt(P_tR_all(tt,:) * (mu_all(tt) - t_roi_decoding') .^ 2);
        end
        
        % compute posterior median
        median_flags_ramp = [false(T_roi_decoding,1),diff(cumsum(P_tR_ramp,2) > .5,1,2) == 1];
        median_flags_non = [false(T_roi_decoding,1),diff(cumsum(P_tR_non,2) > .5,1,2) == 1];
%         median_flags_all = [false(T_roi_decoding,1),diff(cumsum(P_tR_all,2) > .5,1,2) == 1];
        [~,median_idcs_ramp] = max(median_flags_ramp,[],2);
        [~,median_idcs_non] = max(median_flags_non,[],2);
%         [~,median_idcs_all] = max(median_flags_all,[],2);
        med_ramp = t_roi_decoding(median_idcs_ramp);
        med_non = t_roi_decoding(median_idcs_non);
%         med_all = t_roi_decoding(median_idcs_all);
        
        % compute posterior IQR
        q25_flags_ramp = [false(T_roi_decoding,1),diff(cumsum(P_tR_ramp,2) > .25,1,2) == 1];
        q75_flags_ramp = [false(T_roi_decoding,1),diff(cumsum(P_tR_ramp,2) > .75,1,2) == 1];
        q25_flags_non = [false(T_roi_decoding,1),diff(cumsum(P_tR_non,2) > .25,1,2) == 1];
        q75_flags_non = [false(T_roi_decoding,1),diff(cumsum(P_tR_non,2) > .75,1,2) == 1];
%         q25_flags_all = [false(T_roi_decoding,1),diff(cumsum(P_tR_all,2) > .25,1,2) == 1];
%         q75_flags_all = [false(T_roi_decoding,1),diff(cumsum(P_tR_all,2) > .75,1,2) == 1];
        [~,q25_idcs_ramp] = max(q25_flags_ramp,[],2);
        [~,q75_idcs_ramp] = max(q75_flags_ramp,[],2);
        [~,q25_idcs_non] = max(q25_flags_non,[],2);
        [~,q75_idcs_non] = max(q75_flags_non,[],2);
%         [~,q25_idcs_all] = max(q25_flags_all,[],2);
%         [~,q75_idcs_all] = max(q75_flags_all,[],2);
        iqr_ramp = t_roi_decoding(q75_idcs_ramp) - t_roi_decoding(q25_idcs_ramp);
        iqr_non = t_roi_decoding(q75_idcs_non) - t_roi_decoding(q25_idcs_non);
%         iqr_all = t_roi_decoding(q75_idcs_all) - t_roi_decoding(q25_idcs_all);
        
        % compute MAP
        [~,mode_idcs_ramp] = max(P_tR_ramp,[],2);
        [~,mode_idcs_non] = max(P_tR_non,[],2);
%         [~,mode_idcs_all] = max(P_tR_all,[],2);
        map_ramp = t_roi_decoding(mode_idcs_ramp);
        map_non = t_roi_decoding(mode_idcs_non);
%         map_all = t_roi_decoding(mode_idcs_all);
        
        %% store current simulation
        P_TR_RAMP(:,:,ss,mm) = P_tR_ramp;
        P_TR_NON(:,:,ss,mm) = P_tR_non;
%         P_TR_ALL(:,:,ss,mm) = P_tR_all;
        MAP_ramp(:,ss,mm) = map_ramp;
        MAP_non(:,ss,mm) = map_non;
%         MAP_all(:,ss,mm) = map_all;
        MU_ramp(:,ss,mm) = mu_ramp;
        MU_non(:,ss,mm) = mu_non;
%         MU_all(:,ss,mm) = mu_all;
        MED_ramp(:,ss,mm) = med_ramp;
        MED_non(:,ss,mm) = med_non;
%         MED_all(:,ss,mm) = med_all;
        SD_ramp(:,ss,mm) = sd_ramp;
        SD_non(:,ss,mm) = sd_non;
%         SD_all(:,ss,mm) = sd_all;
        IQR_ramp(:,ss,mm) = iqr_ramp;
        IQR_non(:,ss,mm) = iqr_non;
%         IQR_all(:,ss,mm) = iqr_all;
        MUS(:,ss,mm) = mus(selected_idcs);
        SIGMAS(:,ss,mm) = sigmas(selected_idcs);
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
            train_flags = ismember(1:K,randperm(K,K/2));
            rho = corrcoef(...
                nanmean(R(:,nn,train_flags),3),...
                nanmean(R(:,nn,~train_flags),3));
            stereotypy(nn) = rho(1,2);
        end
        
        % store metrics
        TUNING(:,ss,mm) = tuning(selected_idcs);
        FR_RANGE(:,ss,mm) = fr_range(selected_idcs);
        STEREOTYPY(:,ss,mm) = stereotypy(selected_idcs);
        
        %% skip plotting single-run stuff if it's not the last run
        if ss < S
            continue;
        end
        
        %% tiling
        
        % figure initialization
        figure(figopt,...
            'windowstyle','docked',...
            'name',sprintf('tiling_%i',mm));
        
        % axes initialization
        n_rows = 2;
        n_cols = 3;
        n_sps = n_rows * n_cols;
        sps = gobjects(n_sps,1);
        for ii = 1 : n_sps
            sps(ii) = subplot(n_rows,n_cols,ii);
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
        set(sps(1+[0,n_cols]),...
            'ylim',[1,N_clus*2],...
            'ytick',[1,N_clus*2]);
        set(sps([2;3]+[0,n_cols]),...
            'ylim',[1,N_clus],...
            'ytick',[1,N_clus])
        title(sps(1),'All');
        title(sps(2),'Ramps');
        title(sps(3),'Non-ramps');
        
        % normalization
        mus = nanmean(r);
        sigs = nanstd(r);
        z = (r - mus) ./ sigs;

        % color limits
        r_clim = [min(r,[],'all'),max(r,[],'all')];
        z_clim = [-1,1] .* max(abs(z),[],'all');
        
        % plot rate rasters
        imagesc(sps(1),t,[1,N_clus*2],r(:,selected_idcs)',r_clim);
        imagesc(sps(2),t,[1,N_clus],...
            r(:,all_ramp_idcs(selected_ramp_idcs))',r_clim);
        imagesc(sps(3),t,[1,N_clus],...
            r(:,all_non_idcs(selected_non_idcs))',r_clim);
        
        % plot z-score rasters
        imagesc(sps(4),t,[1,N_clus*2],z(:,selected_idcs)',z_clim);
        imagesc(sps(5),t,[1,N_clus],...
            z(:,all_ramp_idcs(selected_ramp_idcs))',z_clim);
        imagesc(sps(6),t,[1,N_clus],...
            z(:,all_non_idcs(selected_non_idcs))',z_clim);

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
                xlabel(sps(rr,cc),'Real time (ms)');
                ylabel(sps(rr,cc),'Decoded time (ms)');
            end
        end
        set(sps,...
            'xlim',roi_decoding,...
            'ylim',roi_decoding,...
            'xdir','normal',...
            'ydir','normal',...
            'nextplot','add',...
            'plotboxaspectratio',[1,1,1]);
        linkaxes(sps);
        
        % iterate through conditions
        clim = quantile([P_tR_ramp(:);P_tR_non(:)],[.001,.999]);
        
        % ramping posteriors
        title(sps(1),sprintf('Ramping neurons (%i/%i)',sum(ramp_flags),N));
        imagesc(sps(1),roi_decoding,roi_decoding,P_tR_ramp',clim);
        plot(sps(1),xlim(sps(1)),ylim(sps(1)),'-k');
        plot(sps(1),xlim(sps(1)),ylim(sps(1)),'--w');
        
        % non-ramping posteriors
        title(sps(2),sprintf('Non-ramping neurons (%i/%i)',sum(~ramp_flags),N));
        imagesc(sps(2),roi_decoding,roi_decoding,P_tR_non',clim);
        plot(sps(2),xlim(sps(2)),ylim(sps(2)),'-k');
        plot(sps(2),xlim(sps(2)),ylim(sps(2)),'--w');
        
        % all posteriors
%         title(sps(3),sprintf('All neurons (%i/%i)',N_clus,N));
%         imagesc(sps(3),roi_decoding,roi_decoding,P_tR_all',clim);
%         plot(sps(3),xlim(sps(3)),ylim(sps(3)),'-k');
%         plot(sps(3),xlim(sps(3)),ylim(sps(3)),'--w');
    end
end

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
bounds.eta = [0,range(roi)*3/4];

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
    'xtick',unique([roi_onset,roi,roi_offset,t_set']));

% within-epoch unimodality assessment of temporal tuning
rois = [roi_onset; roi; roi_offset];
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
bounds.tuning = roi;
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
    'xtick',unique([roi_onset,roi,roi_offset,t_set']));

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

%  model_labels = {...
%     '$\eta=\eta_0\newline\gamma=\gamma_0$';...
%     '$\eta\sim\newline\gamma=\gamma_0$';...
%     '$\eta=\eta_0\newline\gamma\sim$';...
%     '$\eta\sim\newline\gamma\sim$';...
%     };

% axes initialization
xxmax = M + 2;
xxtick = unique((1:xxmax)+[-1;0;1]*.05*xxmax);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:xxmax)) = {''};
xxticklabel(ismember(xxtick,1:M)) = model_labels;
xxticklabel(ismember(xxtick,M+[1,2])) = {'S1';'S2'};
axes(axesopt.default,...
    'plotboxaspectratio',[2.5,1,1],...
    ...'ticklabelinterpreter','latex',...
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
%     all_sims(:,mm) = accuracyfun(pthat_all(:,:,mm),1);
    ramp_avg = avgfun(ramp_sims(:,mm));
    non_avg = avgfun(non_sims(:,mm));
%     all_avg = avgfun(all_sims(:,mm));
    ramp_err = errfun(ramp_sims(:,mm));
    non_err = errfun(non_sims(:,mm));
%     all_err = errfun(all_sims(:,mm));
    plot(mm+xoffsets,[ramp_avg,non_avg],...,all_avg],...
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
%     errorbar(mm+xoffsets(3),all_avg,...
%         all_err(1),all_err(2),...
%         'color','k',...
%         'marker','o',...
%         'markersize',7.5,...
%         'markeredgecolor','k',...
%         'markerfacecolor','w',...
%         'linewidth',1.5,...
%         'capsize',0);
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

% iterate through models
for mm = 1 : M
    xx = [-1,1] * .5 / 3 + mm;
    yy = [1,1] * max(yylim);
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5,...
        'handlevisibility','off');
    [~,pval] = kstest2(ramp_sims(:,mm),non_sims(:,mm));
    pval = kruskalwallis([ramp_sims(:,mm),non_sims(:,mm)],[],'off');
    pval = pval * M;
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
%     all_sims(:,mm) = nanmean(errhat_all(:,:,mm),1);
    ramp_avg = avgfun(ramp_sims(:,mm));
    non_avg = avgfun(non_sims(:,mm));
%     all_avg = avgfun(all_sims(:,mm));
    ramp_err = errfun(ramp_sims(:,mm));
    non_err = errfun(non_sims(:,mm));
%     all_err = errfun(all_sims(:,mm));
    plot(mm+xoffsets,[ramp_avg,non_avg],...,all_avg],...
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
%     errorbar(mm+xoffsets(3),all_avg,...
%         all_err(1),all_err(2),...
%         'color','k',...
%         'marker','o',...
%         'markersize',7.5,...
%         'markeredgecolor','k',...
%         'markerfacecolor','w',...
%         'linewidth',1.5,...
%         'capsize',0);
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
    pval = kruskalwallis([ramp_sims(:,mm),non_sims(:,mm)],[],'off');
    pval = pval * M;
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
            xlabel(sps(rr,cc),'Time (ms)');
            ylabel(sps(rr,cc),'Decoded time (ms)');
        end
    end
    set(sps,...
        axesopt.default,...
        'xlim',roi_decoding,...
        'ylim',roi_decoding,...
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
    'xlim',roi_decoding,...
    'ylim',roi_decoding,...
    'xtick',unique([roi_decoding,t_set']),...
    'ytick',unique([roi_decoding,t_set']),...
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

% all posteriors
% all_avg = nanmean(pthat_all(:,:,M),2);
% all_err = [1,1] .* nanmean(errhat_all(:,:,M),2);
% xpatch = [t_roi_decoding,fliplr(t_roi_decoding)];
% ypatch = [all_avg-all_err(:,1);flipud(all_avg+all_err(:,2))];
% patch(xpatch,ypatch,'w',...
%     'facealpha',1,...
%     'edgecolor','none');

% plot averages
plot(t_roi_decoding,non_avg,...
    'color',([1,1,1]+ramp_clrs(2,:))/2,...
    'linewidth',1.5);
plot(t_roi_decoding,ramp_avg,...
    'color',([1,1,1]+ramp_clrs(1,:))/2,...
    'linewidth',1.5);
% plot(t_roi_decoding,all_avg,...
%     'color',[1,1,1],...
%     'linewidth',1.5);

% inset with example run
axes(axesopt.default,...
    axesopt.inset.nw,...
    'xaxislocation','bottom',...
    'xlim',roi_decoding,...
    'ylim',roi_decoding,...
    'xtick',roi_decoding,...
    'ytick',roi_decoding);

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

%%
figure;
for mm = 1 : M
    subplot(M,1,mm);
    histogram(P_RAMP(:,mm),(0:.025:1.05)-.0125);
    xlim([0,1]);
end
figure(...
    'position',[1.0258e+03 41.8000 1.0224e+03 1.0288e+03]);
fr_edges = linspace(0,30,50);
for mm = 1 : M
    subplot(M,1,mm);
    set(gca,axesopt.default,...
        'plotboxaspectratio',[M,1,1]);
    xlabel('Firing rate (Hz)');
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
end
return;

%% model specification

% figure initialization
fig = figure(figopt,...
    'position',[744,630,415,420],...
    'name','ramps_model_specification');

% model spec
mu_mdl = .5 * range(roi);
gamma_mdl = 5;
sigma_mdl = .15 * range(roi);
x_mdl = generativerate(t_roi,gamma_mdl,mu_mdl,0,sigma_mdl);

% axes initialization
axes(axesopt.default,...
    'xcolor','k',...
    'xlim',roi+[-1,1]*.05*range(roi),...
    'xtick',unique([roi,mu_mdl+[-1,0,1]*sigma_mdl]),...
    'xticklabel',{num2str(roi(1)),'$\mu-\sigma$','$\mu$','$\mu+\sigma$',num2str(roi(2))},...
    'ylim',[0,gamma_mdl]+[-1,1]*.15*range([0,gamma_mdl]),...
    'ytick',[0,gamma_mdl],...
    'yticklabel',{'$b$','$\gamma$'},...
    'ticklabelinterpreter','latex',...
    'clipping','off',...
    'plotboxaspectratio',[3,1,1]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');

% model model illustration
plot(t_roi,x_mdl,...
    'color','k',...
    'linewidth',1.5);

% annotate model spcification
text(.5,1.25,...
    ['$\rm{FR}=\gamma\cdot\mathcal{N}',...
    '(\nu\sim\mathcal{N}(\mu,\eta),\sigma)$'],...
    'color','k',...
    'fontsize',12,...
    'interpreter','latex',...
    'horizontalalignment','center',...
    'units','normalized');
text(.5,1.05,...
    ['$\lambda_{n,k}(t)\propto\gamma_{n}\cdot\mathcal{N}',...
    '(\nu_{k}\sim\mathcal{N}(\mu_{n},\eta_{n}),\sigma_{n})$'],...
    'color','k',...
    'fontsize',12,...
    'interpreter','latex',...
    'horizontalalignment','center',...
    'units','normalized');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% example neurons with opposite gammas & etas
K_eg = 10;

% preallocation
N_eg = 2;
X_eg = nan(T,N_eg,K_eg);
R_eg = nan(T,N_eg,K_eg);

% parameter choices
mus_eg = [...
    [1,1].*.75;...
    [1,1].*.25] * tf;
etas_eg = [...
    [0,0];...
    [0,1]*.25] * tf;
gammas_eg = [...
    [4,2];...
    [3,3]];
sigma_eg = .15 * (tf - ti);
n_eg = size(mus_eg,1);

% iterate through examples
for ee = 1 : n_eg
    
    % figure initialization
    fig = figure(figopt,...
        'position',[500+(ee-1)*300,350,200,275],...
        'name',sprintf('ramps_simulation_eg%i',ee));
    
    % axes settings
    yspacer = mean(gammas_eg(ee,:));
    axes(axesopt.default,...
        'xlim',[ti,tf],...
        'xtick',[ti,tf],...
        'ytick',[1,K_eg*yspacer],...
        'ytick',(1:K_eg)*yspacer,...
        'yticklabel',num2cell(1:K_eg),...
        'ylim',[1,K_eg+1]*yspacer+[-1,0]*.05*K_eg*yspacer,...
        'ylimspec','tight',...
        'ycolor','none',...
        'clipping','off',...
        'ticklength',axesopt.default.ticklength*2,...
        'plotboxaspectratio',[1,1,1]);
    xlabel('Time (ms)');
    ylabel('Trials');
    
    % iterate through neurons
    for nn = 1 : N_eg
        
        % iterate through trials
        for kk = 1 : K_eg
            X_eg(:,nn,kk) = generativerate(...
                t,gammas_eg(ee,nn),mus_eg(ee,nn),etas_eg(ee,nn),sigma_eg);
            [n,ts] = poissonprocess(X_eg(:,nn,kk),(tf-ti)/t_units);
            spk_times = ts * t_units + ti;
            spk_counts = histcounts(spk_times,...
                'binedges',[t,t(end)+1]);
            R_eg(:,nn,kk) = movsum(spk_counts,dt) / (dt/t_units);
        end
    end
    
    % compute cross-trial mean
    x_eg = nanmean(X_eg,3);
    
    % iterate through trials
    for kk = K_eg : -1 : 1
        
        % iterate through neurons
        for nn = 1 : N_eg
            
            % plot single trial
            plot(t,X_eg(:,nn,kk)+kk*yspacer,...
                'color','w',...ramp_clrs(nn,:),...
                'linewidth',1.5);
            xpatch = [fliplr(t),t];
            ypatch = [zeros(T,1);X_eg(:,nn,kk)] + kk * yspacer;
            patch(xpatch,ypatch,ramp_clrs(nn,:),...
                'facealpha',1,...
                'edgecolor','none',...
                'linewidth',1);
        end
    end
    
    % parameter annotation
    if ee == 1
        text(.5,1.3,'$\mu_1=\mu_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.2,'$\gamma_1>\gamma_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.1,'$\eta_1=\eta_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
    else
        text(.5,1.3,'$\mu_1=\mu_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.2,'$\gamma_1=\gamma_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.1,'$\eta_1>\eta_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
    end
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end

%% firing rate function
function x = generativerate(time,gamma,mu,eta,sigma)
    mu_sample = clamp(normrnd(mu,eta),time(1),time(end));
    x_pdf = normpdf(time,mu_sample,sigma);
    x = gamma * x_pdf ./ max(x_pdf);
end