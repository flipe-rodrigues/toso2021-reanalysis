%% initialization
if ~exist('data','var')
    toso2021_wrapper;
    close all;
end

%% training & test set conditions

% training set conditions
if strcmpi(contrast_str,'t1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'choice')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'correct')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'correct',correct(valid_flags),1,[]);
elseif strcmpi(contrast_str,'choice_correct')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[],...
        'correct',correct(valid_flags),[],[]);
end

% test set conditions
if strcmpi(contrast_str,'t1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),t_set,[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i1')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set,[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'choice')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[0,1],[]);
elseif strcmpi(contrast_str,'correct')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...
        'correct',correct(valid_flags),[0,1],[]);
elseif strcmpi(contrast_str,'choice_correct')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...
        'i2',i2(valid_flags),[],[],...,...
        'choice',choice(valid_flags),[0,1],[],...
        'correct',correct(valid_flags),[0,1],[]);
end
n_conditions = conditions.test.n + conditions.train.n;

% print training & test conditions
fprintf('\nTRAINING CONDITIONS:\n');
conditions.train.values
fprintf('\nTEST CONDITIONS:\n');
conditions.test.values

%% run settings
n_runs = 100;

%% subject selection
subject_flags = ismember(subjects,subject_set);

%% concatenation settings
n_concats_max = 2^7;

%% time settings
roi = [-500,t_set(end)];
roi2plot = [-240,t_set(t2_mode_idx+1)];
roi_n_bins = range(roi) / psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);
roi_xlim = [-240,t_set(t2_mode_idx+1)]; % [-350,t_set(t2_mode_idx+2)];
roi_ylim = [-240,t_set(t2_mode_idx+1)];

%% slice settings
slice_times = unique([roi2plot(1)/2;0;t_set(1:end-2)]);
n_slices = numel(slice_times);
slice_offsets = (n_slices - 1 : -1 : 0) * .05 + .015;
slice_clrs = gray(n_slices);

%% construct spike rate tensor (time X neurons X concatenations)

% data clearance
clear R P_tR;

% preallocation
P_tR = nan(roi_n_bins,roi_n_bins,conditions.test.n,n_runs);
map = nan(roi_n_bins,conditions.test.n,n_runs);

% temporal smoothing kernel
gauss_kernel = gausskernel('sig',50,'binwidth',psthbin);

% iterate through runs
for rr = 1 : n_runs
    
    % preallocation
    R = nan(roi_n_bins,n_neurons,n_conditions);

    % iterate through units
    for nn = 1 : n_neurons
        progressreport(nn,n_neurons,...
            sprintf('sampling concatenations (run %i/%i)',rr,n_runs));
        neuron_flags = data.NeuronNumb == flagged_neurons(nn);
        
        % preallocation
        xval_train_trials = cell(conditions.train.n,conditions.test.n);
        
        % iterate through training conditions
        for kk = 1 : conditions.train.n
            
            % flag trials for the current condition
            feature_flags = false(n_total_trials,conditions.train.features.n);
            for ff = 1 : conditions.train.features.n
                feature_lbl = conditions.train.features.labels{ff};
                feature = eval(feature_lbl);
                feature_flags(:,ff) = ismember(...
                    feature,conditions.train.values.(feature_lbl)(kk,:));
            end
            condition_flags = all(feature_flags,2);
            
            % trial selection
            trial_flags = ...
                valid_flags & ...
                subject_flags & ...
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_rates = data.SDF(trial_flags,:);
            
            % S2-onset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi;
            alignment_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + t2(trial_flags);
            chunk_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + roi(2);
            aligned_spkrates = spike_rates';
%             aligned_spkrates = ...
%                 conv2(gauss_kernel.pdf,1,aligned_spkrates,'same');
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [roi_n_bins,n_flagged_trials])';
            
            % detect any test conditions that overlap with the current training condition
            feature_flags = false(conditions.test.n,1);
            for ff = 1 : conditions.test.features.n
                feature_lbl = conditions.test.features.labels{ff};
                feature_flags(:,ff) = any(ismember(...
                    conditions.test.values.(feature_lbl),...
                    conditions.train.values.(feature_lbl)(kk,:)),2);
            end
            xval_condition_flags = all(feature_flags,2);
            
            % handle cross-validation for all detected test conditions
            if any(xval_condition_flags)
                xval_condition_idcs = find(xval_condition_flags)';
                
                % iterate through detected test conditions
                for ii = xval_condition_idcs
                    
                    % detect overlap between training & test trials
                    feature_flags = false(n_flagged_trials,1);
                    for ff = 1 : conditions.test.features.n
                        feature_lbl = conditions.test.features.labels{ff};
                        feature = eval(feature_lbl);
                        feature_flags(:,ff) = any(ismember(...
                            feature(flagged_trials),...
                            conditions.test.values.(feature_lbl)(ii,:)),2);
                    end
                    xval_trial_flags = all(feature_flags,2);
                    
                    % split the conflicting trials into training & test subsets
                    xval_trials = flagged_trials(xval_trial_flags);
                    n_xval_trials = numel(xval_trials);
                    n_train_trials = round(n_xval_trials * 2 / 3);
                    xval_train_idcs = randperm(n_xval_trials,n_train_trials);
                    xval_train_trials{kk,ii} = xval_trials(xval_train_idcs);
                end
                
                % concatenate sub-sampled training sets across test conditions
                train_idcs = find(ismember(...
                    flagged_trials,vertcat(xval_train_trials{kk,:})));
            else
                
                % train using all trials in the remaining conditions
                train_idcs = 1 : n_flagged_trials;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(train_idcs,n_concats_max,true);
            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
            nan_flags = isnan(r_mu);

%             mdl_spline = fit(...
%                 roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
%                 'smoothingparam',1e-6);
%             r_spline = max(1e-3,mdl_spline(roi_time));
%             r_spline(nan_flags) = nan;
%             
%             t_mat = repmat(roi_time,n_concats_max,1)';
%             r_mat = r';
%             r_vec = r_mat(~isnan(r_mat));
%             t_vec = t_mat(~isnan(r_mat));
%             mdl_poly = fit(t_vec,r_vec,'poly9');
%             r_poly = max(realmin,mdl_poly(roi_time));
%             r_poly(nan_flags) = nan;
            
%             r_gauss = nanconv2(r_mu,1,gauss_kernel.pdf);
%             r_gauss(nan_flags) = nan;
            
            r_gauss = conv2(1,gauss_kernel.pdf,r_mu,'valid');
            r_gauss = padarray(r_gauss,[0,floor(gauss_kernel.nbins/2)],nan,'pre');
            r_gauss = padarray(r_gauss,[0,ceil(gauss_kernel.nbins/2)-1],nan,'post');
            r_gauss(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
% %             plot(roi_time,r_spline);
%             plot(roi_time,r_gauss,'k');
            
            R(:,nn,kk) = r_gauss;
        end
        
        % iterate through conditions
        for kk = 1 : conditions.test.n
            
            % flag trials for the current condition
            feature_flags = false(n_total_trials,conditions.test.features.n);
            for ff = 1 : conditions.test.features.n
                feature_lbl = conditions.test.features.labels{ff};
                feature = eval(feature_lbl);
                feature_flags(:,ff) = ismember(...
                    feature,conditions.test.values.(feature_lbl)(kk,:));
            end
            condition_flags = all(feature_flags,2);
            
            % trial selection
            trial_flags = ...
                valid_flags & ...
                subject_flags & ...
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_rates = data.SDF(trial_flags,:);
            
            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi;
            alignment_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + t2(trial_flags);
            chunk_flags = ...
                padded_time >= alignment + roi(1) & ...
                padded_time < alignment + roi(2);
            aligned_spkrates = spike_rates';
%             aligned_spkrates = ...
%                 conv2(gauss_kernel.pdf,1,aligned_spkrates,'same');
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [roi_n_bins,n_flagged_trials])';
            
            % handle cross-validation
            test_idcs = find(~ismember(...
                flagged_trials,vertcat(xval_train_trials{:,kk})));
            if isempty(test_idcs)
                continue;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(test_idcs,n_concats_max,true);
            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
            nan_flags = isnan(r_mu);
            
%             mdl_spline = fit(...
%                 roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
%                 'smoothingparam',1e-6);
%             r_spline = max(1e-3,mdl_spline(roi_time));
%             r_spline(nan_flags) = nan;
            
%             t_mat = repmat(roi_time,n_concats_max,1)';
%             r_mat = r';
%             r_vec = r_mat(~isnan(r_mat));
%             t_vec = t_mat(~isnan(r_mat));
%             mdl_poly = fit(t_vec,r_vec,'poly9');
%             r_poly = max(realmin,mdl_poly(roi_time));
%             r_poly(nan_flags) = nan;
            
%             r_gauss = nanconv2(r_mu,1,gauss_kernel.pdf);
%             r_gauss(nan_flags) = nan;
            
            r_gauss = conv2(1,gauss_kernel.pdf,r_mu,'valid');
            r_gauss = padarray(r_gauss,[0,floor(gauss_kernel.nbins/2)],nan,'pre');
            r_gauss = padarray(r_gauss,[0,ceil(gauss_kernel.nbins/2)-1],nan,'post');
            r_gauss(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
% %             plot(roi_time,r_spline);
%             plot(roi_time,r_gauss,'--k');
            
            R(:,nn,kk+conditions.train.n) = r_gauss;
        end
    end
    
    %%
%     figure;
%     for nn2 = 1 : n_neurons
%         plot(roi_time,R(:,nn2,1))
%         title(sprintf('%i',nn2));
%         pause(.01);
%         if any(isnan(squeeze(R(:,nn2,1))))
%             nn2
%         end
%     end
    
    %% naive bayes decoder
    nbdopt = struct();
    nbdopt.n_xpoints = 100;
    nbdopt.time = roi_time;
    nbdopt.train.trial_idcs = 1 : conditions.train.n;
    nbdopt.train.n_trials = numel(nbdopt.train.trial_idcs);
    nbdopt.test.trial_idcs = (1 : conditions.test.n) + conditions.train.n;
    nbdopt.test.n_trials = numel(nbdopt.test.trial_idcs);
    nbdopt.assumepoissonmdl = true;
    nbdopt.verbose = true;
    
    tic
    [P_tR(:,:,:,rr),P_Rt,pthat,neurons] = naivebayestimedecoder(R,nbdopt);
    map(:,:,rr) = pthat.mode;
    toc
end

%% plot likelihoods
if false
    figure;
    set(gca,...
        axesopt.default,...,...
        'xlim',roi,...
        'xtick',sort([0;t_set]));
    xlabel('Time since T_2 (ms)');
    ylabel('Firing rate (Hz)');
    
    % iterate through units
    for nn = 1 : n_neurons
        cla;
        r_bounds = neurons(nn).x_bounds;
        r_bw = neurons(nn).x_bw;
        if range(r_bounds) == 0
            continue;
        end
        ylim(r_bounds);
        title(sprintf('neuron: %i, bw: %.2f',nn,r_bw));
        p_Rt = squeeze(P_Rt(:,nn,:));
        nan_flags = isnan(p_Rt);
        if sum(abs(diff(any(nan_flags,2)))) > 1
            fprintf('\tcheck neuron %i!\n',nn);
        end
        p_Rt(nan_flags) = max(p_Rt(:));
        imagesc(xlim,r_bounds,p_Rt');
        for ii = 1 : n_t
            plot([1,1]*t_set(ii),ylim,'--w');
        end
        drawnow;
        pause(.1);
    end
end

%% choice of average function
avgfun = @(x,d)nanmean(x,d);
errfun = @(x,d)nanstd(x,0,d);

%% plot contrast-split posterior averages
figure(...
    'name','contrast-split posterior averages',...
    'numbertitle','off',...
    'windowstyle','docked');
sps = gobjects(n_contrasts,1);
for ii = 1 : n_contrasts
    sps(ii) = subplot(1,n_contrasts,ii);
    xlabel(sps(ii),'Time since S_2 onset (ms)');
    ylabel(sps(ii),'Decoded time since S_2 onset (ms)');
end
set(sps,...
    axesopt.default,...
    'xlim',roi_xlim,...
    'xtick',unique([roi';roi_xlim';0;t_set]),...
    'ylim',roi_ylim,...
    'ytick',unique([roi';roi_ylim';0;t_set]));
linkaxes(sps);

clims = quantile(P_tR,[0,.999],'all')';

% iterate through contrast conditions
for ii = 1 : n_contrasts
    title(sps(ii),sprintf('%s = %.0f mm/s',...
        contrast_lbl,contrast_set(ii)),...
        'fontsize',10);
    p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
    imagesc(sps(ii),roi,roi,p_cond',clims);
    plot(sps(ii),xlim,ylim,'--w');
end

%% plot superimposed contrast-split posterior averages
fig = figure(...
    figopt,...
    'name',sprintf('superimposed_posterior_averages_%s',contrast_str),...
    'numbertitle','off');
axes(...
    axesopt.default,...
    'xlim',roi_xlim,...+[-1,1]*.05*range(roi2plot),...
    'xtick',unique([roi';roi_xlim';0;t_set]),...
    'ylim',roi_ylim,...+[-1,1]*.05*range(roi2plot),...
    'ytick',unique([roi';roi_ylim';0;t_set]),...
    'xticklabelrotation',0,...
    'yticklabelrotation',0);
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% convert from tensor to rgb
P_tR_avg = squeeze(avgfun(P_tR,4));
% P_tR_avg(roi_time < roi_time(1) | roi_time > roi2plot(2),:,:) = nan;
% P_tR_avg(:,roi_time < roi_time(1) | roi_time > roi2plot(2),:) = nan;
% P_tR_avg(isnan(P_tR_avg)) = 0;
% P_tR_avg = log(P_tR_avg);

% P_tR_avg = min(P_tR_avg,quantile(P_tR_avg,.9975,[1,2]));

P = tensor2rgb(permute(P_tR_avg,[2,1,3]),contrast_clrs);
imagesc(roi,roi,P);

% zero lines
plot(xlim,[1,1]*0,':k');

% plot identity line
plot([roi2plot(1),t_set(end)],[roi2plot(1),t_set(end)],':k');

% iterate through slices
for ii = 1 : n_slices
    
    % plot slice handles
    plot([1,1]*slice_times(ii),ylim,':k');
    plot(slice_times(ii),max(ylim),...
        'marker','v',...
        'markersize',7.5,...
        'markerfacecolor',slice_clrs(ii,:),...
        'markeredgecolor','k',...
        'linewidth',1.5);
end

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
    [contrast_clrs(1,:);[1,1,1];contrast_clrs(end,:)],2^8));
ylabel('P(t|R)',...
    'verticalalignment','middle',...
    'rotation',-90);

% colorbar settings
clrbar_width = .05;

% iterate through contrasts
for ii = 1 : n_contrasts
    
    % patch pseudo-colorbar
    xpatch = (1 - clrbar_width * n_contrasts) + ...
        clrbar_width * ((ii - 1) + [0,1,1,0]);
    ypatch = [0,0,1,1];
    patch(xpatch,ypatch,contrast_clrs(ii,:),...
        'edgecolor','none',...
        'linewidth',1.5);
end

% inset with extreme posterior subtractions
% axes(...
%     axesopt.inset.se,...
%     axesopt.default,...
% 	'xlim',roi_xlim,...roi2plot,...+[-1,1]*.05*range(roi2plot),...
%     'xtick',unique([roi';roi_xlim';0;t_set]),...
%     'ylim',roi_ylim,...+[-1,1]*.05*range(roi2plot),...
%     'ytick',unique([roi';roi_ylim';0;t_set]),...
%     'box','on',...
%     'colormap',colorlerp(...
%     [contrast_clrs(1,:);[1,1,1];contrast_clrs(end,:)],2^8));
% % title('$p(t|\textbf{r})-p(t|\textbf{r})$',...
% %     'interpreter','latex');
% title('P(t|R) - P(t|R)');
% ylabel('\DeltaP(t|R)',...
%     'rotation',-90,...
%     'verticalalignment','bottom');
% 
% % posterior subtraction
% p_contrast_min = avgfun(P_tR(:,:,1,:),4);
% p_contrast_max = avgfun(P_tR(:,:,end,:),4);
% p_contrast_min = p_contrast_min - min(p_contrast_min,[],[1,2]);
% p_contrast_max = p_contrast_max - min(p_contrast_max,[],[1,2]);
% p_diff = p_contrast_max - p_contrast_min;
% imagesc(roi,roi,p_diff',[-1,1] * n_t / n_tbins * 5);
% 
% % plot identity line
% plot([roi2plot(1),t_set(end)],[roi2plot(1),t_set(end)],':k');
% 
% % zero lines
% plot(xlim,[1,1]*0,':k');
% plot([1,1]*0,ylim,':k');

% inset with contrast-split MAPs
axes(...
    axesopt.inset.se,...
    axesopt.default,...
    'box','off',...
    'color','w',...
    'xlim',[0,roi_xlim(2)],...
    'xtick',unique([roi';roi2plot';0;t_set]),...
    'ylim',[0,roi_ylim(2)],...
    'ytick',unique([roi';roi2plot';0;t_set]));
% title('M.A.P.');

% iterate through contrast conditions
for ii = 1 : n_contrasts
    map_avg = squeeze(avgfun(map(:,ii,:),3));
    map_err = squeeze(errfun(map(:,ii,:),3));
    errorpatch(roi_time,map_avg,map_err,contrast_clrs(ii,:),...
        'facealpha',.25);
    plot(roi_time,map_avg,...
        'color',contrast_clrs(ii,:),...
        'linewidth',1.5);
end

% % zero lines
% plot([1,1]*0,ylim,':k');
% plot(xlim,[1,1]*0,':k');
% 
% % plot identity line
% plot([roi2plot(1),t_set(end)],[roi2plot(1),t_set(end)],':k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot superimposed contrast-split MAP averages
fig = figure(...
    figopt,...
    'name',sprintf('superimposed_posterior_maps_%s',contrast_str),...
    'numbertitle','off');
axes(...
    axesopt.default,...
    'xlim',roi_xlim,...+[-1,1]*.05*range(roi2plot),...
    'xtick',unique([roi';roi_xlim';0;t_set]),...
    'ylim',roi_ylim,...+[-1,1]*.05*range(roi2plot),...
    'ytick',unique([roi';roi_ylim';0;t_set]),...
    'xticklabelrotation',0,...
    'yticklabelrotation',0);
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% plot identity line
plot([roi2plot(1),t_set(end)],[roi2plot(1),t_set(end)],':k');

% zero lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*0,':k');

% iterate through contrast conditions
for ii = 1 : n_contrasts
    map_avg = squeeze(avgfun(map(:,ii,:),3));
    map_err = squeeze(errfun(map(:,ii,:),3));
    errorpatch(roi_time,map_avg,map_err,contrast_clrs(ii,:),...
        'facealpha',.25);
    plot(roi_time,map_avg,...
        'color',contrast_clrs(ii,:),...
        'linewidth',1.5);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot slices through condition-split posterior averages

% figure initialization
fig = figure(...
    'color','w',...
    'name',sprintf('posterior_slices_%s',contrast_str),...
    'numbertitle','off');

% axes initialization
axes(...
    axesopt.default,...
    'xlim',roi2plot,...
    'xtick',unique([roi';roi2plot';0;t_set]),...
    'ytick',0,...
    'layer','bottom');

% axes labels
xlabel('Decoded time since S_2 onset (ms)');
ylabel('P(t|R)');

% time seletion
time_flags = ...
    roi_time >= min(xlim) & ...
    roi_time <= max(xlim);

% iterate through contrasts
[~,cond_idcs] = sort(...
    abs(contrast_set-contrast_set(contrast_mode_idx)),'descend');
for ii = cond_idcs'
    p_cond = avgfun(P_tR(:,:,ii,:),4);

    % iterate through slices
    for jj = 1 : n_slices
        slice_idx = find(roi_time >= slice_times(jj),1);
        
        % plot posterior slice
        plot(roi_time(time_flags),p_cond(slice_idx,time_flags)+slice_offsets(jj),...
            'color',contrast_clrs(ii,:),...
            'linewidth',1.5);
    end
end

% iterate through slices
for ii = 1 : n_slices

    % plot real time
    slice_idx = find(t >= slice_times(ii),1);
    plot(min(xlim),slice_offsets(ii),...
        'marker','>',...
        'markersize',7.5,...
        'markerfacecolor',slice_clrs(ii,:),...
        'markeredgecolor','k',...
        'linewidth',1.5);
end

% zero lines
plot([1,1]*0,[0,max(ylim)]*1.5,':k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end