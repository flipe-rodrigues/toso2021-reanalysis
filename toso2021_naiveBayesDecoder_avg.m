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
        'i2',i2(valid_flags),[],[],...t_set(t2_mode_idx+[-1,0,1]),...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.train = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),i_set(i1_mode_idx),[],...
        't2',t2(valid_flags),[],[],...t_set(t2_mode_idx+[-1,0,1]),...
        'i2',i2(valid_flags),i_set(i2_mode_idx),[],...
        'choice',choice(valid_flags),[],[]);
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
        't2',t2(valid_flags),[],[],...t_set([1,2,end-1,end]),...
        'i2',i2(valid_flags),[],[],...
        'choice',choice(valid_flags),[],[]);
elseif strcmpi(contrast_str,'i2')
    conditions.test = intersectconditions(...
        't1',t1(valid_flags),[],[],...
        'i1',i1(valid_flags),[],[],...
        't2',t2(valid_flags),[],[],...t_set([1,2,end-1,end]),...
        'i2',i2(valid_flags),i_set,[],...
        'choice',choice(valid_flags),[],[]);
end

n_conditions = conditions.test.n + conditions.train.n;

% print training & test conditions
fprintf('\nTRAINING CONDITIONS:\n');
conditions.train.values
fprintf('\nTEST CONDITIONS:\n');
conditions.test.values

%% run settings
n_runs = 10;

%% concatenation settings
n_concatspercond = 2^7;
n_concats = n_concatspercond * n_conditions;

%% time settings
roi = [-500,t_set(end)];
roi_n_bins = range(roi) / psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

%% construct spike rate tensor (time X neurons X concatenations)

% data clearance
clear concat_tensor P_tR;

% preallocation
P_tR = nan(roi_n_bins,roi_n_bins,conditions.test.n,n_runs);
MAP = nan(roi_n_bins,conditions.test.n,n_runs);
% concat_stimuli = nan(n_concats,n_runs);
% concat_contrasts = nan(n_concats,n_runs);
% concat_choices = nan(n_concats,n_runs);
% concat_evalset = categorical(nan(n_concats,n_runs),[0,1],{'train','test'});

% data type selection
spike_data_field = 'FR';

% iterate through runs
for rr = 1 : n_runs
    
    % preallocation
    R = nan(roi_n_bins,n_neurons,n_conditions);
    %     concat_tensor = nan(roi_n_bins,n_neurons,n_concats);
    %     concat_stimuli = nan(n_concats,1);
    %     concat_contrasts = nan(n_concats,1);
    %     concat_choices = nan(n_concats,1);
    %     concat_evalset = categorical(nan(n_concats,1),[0,1],{'train','test'});
    
    % iterate through units
    for nn = 1 : n_neurons
        progressreport(nn,n_neurons,...
            sprintf('sampling concatenations (run %i/%i)',rr,n_runs));
        neuron_flags = data.NeuronNumb == flagged_neurons(nn);
        
        % preallocation
        xval_train_trials_bab = cell(conditions.train.n,conditions.test.n);
        
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
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_counts = data.(spike_data_field)(trial_flags,:);
            spike_rates = ...
                conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
            
            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi;
            alignment_flags = ...
                valid_time >= alignment + roi(1) & ...
                valid_time < alignment + t2(trial_flags);
            chunk_flags = ...
                valid_time >= alignment + roi(1) & ...
                valid_time < alignment + roi(2);
            aligned_spkrates = spike_rates;
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
                    xval_train_trials_bab{kk,ii} = xval_trials(xval_train_idcs);
                end
                
                % concatenate sub-sampled training sets across test conditions
                train_idcs = find(ismember(...
                    flagged_trials,vertcat(xval_train_trials_bab{kk,:})));
            else
                
                % train using all trials in the remaining conditions
                train_idcs = 1 : n_flagged_trials;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(train_idcs,n_concatspercond,true);
            %             concat_idcs = (1 : n_concatspercond) + ...
            %                 n_concatspercond * (kk - 1);
            %             concat_tensor(:,nn,concat_idcs) = aligned_spkrates(rand_idcs,:)';
            %             concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
            %             concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
            %             concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
            %             concat_evalset(concat_idcs) = 'train';
            
            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
            nan_flags = isnan(r_mu);

            mdl_spline = fit(...
                roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
                'smoothingparam',1e-6);
            r_spline = max(1e-3,mdl_spline(roi_time));
            r_spline(nan_flags) = nan;
            
%             t_mat = repmat(roi_time,n_concatspercond,1)';
%             r_mat = r';
%             r_vec = r_mat(~isnan(r_mat));
%             t_vec = t_mat(~isnan(r_mat));
%             mdl_poly = fit(t_vec,r_vec,'poly9');
%             r_poly = max(1e-3,mdl_poly(roi_time));
%             r_poly(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
%             plot(roi_time,r_spline);
%             plot(roi_time,(r_poly+r_spline)/2);
            
            R(:,nn,kk) = r_spline;
        end
        
        % iterate through conditions
        for kk = 1 : conditions.test.n
            
            % flag trials for the current condition
            feature_flags = false(n_total_trials,conditions.train.features.n);
            for ff = 1 : conditions.train.features.n
                feature_lbl = conditions.test.features.labels{ff};
                feature = eval(feature_lbl);
                feature_flags(:,ff) = ismember(...
                    feature,conditions.test.values.(feature_lbl)(kk,:));
            end
            condition_flags = all(feature_flags,2);
            
            % trial selection
            trial_flags = ...
                valid_flags & ...
                neuron_flags & ...
                condition_flags;
            flagged_trials = find(trial_flags);
            n_flagged_trials = numel(flagged_trials);
            if n_flagged_trials == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            spike_counts = data.(spike_data_field)(trial_flags,:);
            spike_rates = ...
                conv2(1,kernel.pdf,spike_counts,'valid')' / psthbin * 1e3;
            
            % S2-offset-aligned spike rates
            alignment = ...
                pre_init_padding + ...
                pre_s1_delay(trial_flags) + ...
                t1(trial_flags) + ...
                isi;
            alignment_flags = ...
                valid_time >= alignment + roi(1) & ...
                valid_time < alignment + t2(trial_flags);
            chunk_flags = ...
                valid_time >= alignment + roi(1) & ...
                valid_time < alignment + roi(2);
            aligned_spkrates = spike_rates;
            aligned_spkrates(~alignment_flags') = nan;
            aligned_spkrates = reshape(aligned_spkrates(chunk_flags'),...
                [roi_n_bins,n_flagged_trials])';
            
            % handle cross-validation
            test_idcs = find(~ismember(...
                flagged_trials,vertcat(xval_train_trials_bab{:,kk})));
            if isempty(test_idcs)
                continue;
            end
            
            % store tensor & concatenation data
            rand_idcs = randsample(test_idcs,n_concatspercond,true);
            %             concat_idcs = (1 : n_concatspercond) + ...
            %                 n_concatspercond * (kk + conditions.train.n - 1);
            %             concat_tensor(:,nn,concat_idcs) = aligned_spkrates(rand_idcs,:)';
            %             concat_stimuli(concat_idcs) = stimuli(flagged_trials(rand_idcs));
            %             concat_contrasts(concat_idcs) = contrasts(flagged_trials(rand_idcs));
            %             concat_choices(concat_idcs) = choice(flagged_trials(rand_idcs));
            %             concat_evalset(concat_idcs) = 'test';
            
            r = aligned_spkrates(rand_idcs,:);
            r_mu = nanmean(r);
            nan_flags = isnan(r_mu);
            
            mdl_spline = fit(...
                roi_time(~nan_flags)',r_mu(~nan_flags)','smoothingspline',...
                'smoothingparam',1e-6);
            r_spline = max(1e-3,mdl_spline(roi_time));
            r_spline(nan_flags) = nan;
            
%             t_mat = repmat(roi_time,n_concatspercond,1)';
%             r_mat = r';
%             r_vec = r_mat(~isnan(r_mat));
%             t_vec = t_mat(~isnan(r_mat));
%             mdl_poly = fit(t_vec,r_vec,'poly9');
%             r_poly = max(1e-3,mdl_poly(roi_time));
%             r_poly(nan_flags) = nan;
            
%             figure('position',[119.4000 53.8000 560 712.8000]);
%             subplot(3,1,[1,2]);
%             imagesc(roi_time,[],r);
%             subplot(3,1,3); hold on;
%             plot(roi_time,r_mu);
%             plot(roi_time,r_poly);
%             plot(roi_time,r_spline);
            
            R(:,nn,kk+conditions.train.n) = r_spline;
        end
    end
    
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
    MAP(:,:,rr) = pthat.mode;
    toc
end

%% plot likelihoods
figure;
set(gca,...
    axesopt.default,...,...
    'xlim',roi,...
    'xtick',sort([0;t_set]));
xlabel('Time since T_2 (ms)');
ylabel('Firing rate (Hz)');

% iterate through units
for nn = 263 : n_neurons
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

%% choice of average function
avgfun = @(x,d)nanmean(x,d);
errfun = @(x,d)nanstd(x,0,d);

%% plot extreme subtractions of posterior averages

% figure initialization
fig = figure(...
    figopt,...
    'name','posterior subtractions (extreme)',...
    'numbertitle','off');

% axes initialization
axes(...
    axesopt.default,...
    'xlim',[roi(1),t_set(end-2)],...
    'xtick',unique([roi';0;t_set]),...
    'ylim',[roi(1),t_set(end-2)],...
    'ytick',unique([roi';0;t_set]),...
    'colormap',colorlerp(...
    [contrast_clrs(1,:);[1,1,1];contrast_clrs(end,:)],2^8));
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% posterior subtraction
p_contrast_min = avgfun(P_tR(:,:,1,:),4);
p_contrast_max = avgfun(P_tR(:,:,end,:),4);
p_diff = p_contrast_max - p_contrast_min;
imagesc(roi,roi,p_diff',[-1,1] * n_t / n_tbins * 5);

% zero lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*0,':k');

% identity line
plot(xlim,ylim,'--w');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

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
    'xlim',roi,...
    'xtick',unique([roi';0;t_set]),...
    'ylim',roi,...
    'ytick',unique([roi';0;t_set]));
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
    'color','w',...
    'name','superimposed posterior averages',...
    'numbertitle','off');
axes(...
    axesopt.default,...
    'xlim',[roi(1),t_set(end-2)],...
    'xtick',unique([roi';0;t_set]),...
    'ylim',[roi(1),t_set(end-2)],...
    'ytick',unique([roi';0;t_set]));
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% color limits
clims = quantile(avgfun(P_tR,4),[0,1],'all')';

% zero lines
plot([1,1]*0,ylim,':k');
plot(xlim,[1,1]*0,':k');

% iterate through contrast conditions
for ii = 1 : n_contrasts
    p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
    p_cond(p_cond < clims(1)) = clims(1);
    p_cond(p_cond > clims(2)) = clims(2);
    p_cond(isnan(p_cond)) = max(p_cond(:));
    p_patch = mat2patch(p_cond,roi,roi);
    patch(p_patch,...
        'facevertexalphadata',p_patch.facevertexcdata,...
        'edgecolor','none',...
        'facecolor',contrast_clrs(ii,:),...
        'facealpha','flat');
end

% plot identity line
plot(xlim,ylim,'--k');

% save figure
% if want2save
%     svg_file = fullfile(panel_path,[fig.Name,'.svg']);
%     print(fig,svg_file,'-dsvg','-painters');
% end

%% plot slices through condition-split posterior averages

% figure initialization
fig = figure(...
    'name','slices through condition-split posterior averages',...
    'numbertitle','off',...
    'position',[769.8000 41.8000 766.4000 740.8000],...
    'color','w');

% slice settings
n_slices = 7;
slices = linspace(roi(1),t_set(end-2),n_slices);
% slices = [0,50,100,200,400];
% n_slices = numel(slices);

% axes initialization
sps = gobjects(n_slices,1);
for ii = 1 : n_slices
    sps(ii) = subplot(n_slices,1,ii);
    xlabel(sps(ii),'Decoded time since S_2 onset (ms)');
    ylabel(sps(ii),'P(t|R)');
end
set(sps,...
    'xlim',[roi(1),t_set(end-2)],...
    'xtick',unique([roi';0;t_set]),...
    'ylimspec','tight',...
    'xdir','normal',...
    'ydir','normal',...
    'layer','top',...
    'clipping','on',...
    'fontsize',12,...
    'linewidth',2,...
    'tickdir','out',...
    'nextplot','add',...
    'ticklength',[1,1]*.025,...
    'nextplot','add',...
    'plotboxaspectratio',[4,1,1]);
set(sps,...
    'xlim',xlim(sps(1)) + [-1,1] * .05 * range(xlim(sps(1))));
set(sps(1:end-1),...
    'xcolor','none');

% iterate through contrast conditions
[~,cond_idcs] = sort(...
    abs(contrast_set-contrast_set(contrast_mode_idx)),'descend');
for ii = cond_idcs'
    p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
    
    % iterate through slices
    for jj = 1 : n_slices
        slice_idx = find(roi_time >= slices(jj),1);
        p_slice = p_cond(slice_idx,:);
        p_slice = p_slice / nansum(p_slice);
        
        % plot posterior slice
        plot(sps(jj),roi_time,p_slice,...
            'color',contrast_clrs(ii,:),...
            'linewidth',1.5);
    end
end

% linkage
linkaxes(sps);

% iterate through slices
for ii = 1 : n_slices
    
    % update axes
    set(sps(ii),...
        'ylim',[0,max(ylim(sps(ii)))] + [0,1] * .1 * max(ylim(sps(ii))),...
        'ytick',0,...
        'yticklabel','0');
    
    % plot real time
    slice_idx = find(t >= slices(ii),1);
    plot(sps(ii),t(slice_idx),max(ylim(sps(ii))),...
        'marker','v',...
        'markersize',8.5,...
        'markerfacecolor','k',...
        'markeredgecolor','none');
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot superimposed contrast-split MAPs

% figure initialization
fig = figure(...
    figopt,...
    'name','posterior subtractions (extreme)',...
    'numbertitle','off');

% axes initialization
axes(...
    axesopt.default,...
    'xlim',roi,...
    'xtick',unique([roi';0;t_set]),...
    'ylim',roi,...
    'ytick',unique([roi';0;t_set]));
xlabel('Time since S_2 onset (ms)');
ylabel('Decoded time since S_2 onset (ms)');

% iterate through contrast conditions
for ii = 1 : n_contrasts
    %     p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
    %     [~,mode_idcs] = max(p_cond,[],2);
    %     map_cond = roi_time(mode_idcs);
    map_avg = squeeze(avgfun(MAP(:,ii,:),3));
    map_err = squeeze(errfun(MAP(:,ii,:),3));
    errorpatch(roi_time,map_avg,map_err,contrast_clrs(ii,:),...
        'facealpha',.25);
    plot(roi_time,map_avg,...
        'color',contrast_clrs(ii,:),...
        'linewidth',1.5);
end

% plot identity line
plot(xlim,ylim,'--k');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%%
%
% % figure initialization
% fig = figure(...
%     figopt,...
%     'name','posterior subtractions (extreme)',...
%     'numbertitle','off');
%
% % axes initialization
% axes(...
%     axesopt.default,...
%     'xlim',roi,...
%     'xtick',unique([roi';0;t_set]),...
%     'ylim',roi,...
%     'ytick',unique([roi';0;t_set]));
% xlabel('Time since S_2 onset (ms)');
% ylabel('Decoded time since S_2 onset (ms)');
%
% p_ref = squeeze(avgfun(P_tR(:,:,contrast_mode_idx,:),4));
% [~,mode_idcs] = max(p_ref,[],2);
% map_ref = roi_time(mode_idcs);
%
% % iterate through contrast conditions
% for ii = 1 : n_contrasts
%     p_cond = squeeze(avgfun(P_tR(:,:,ii,:),4));
%     [~,mode_idcs] = max(p_cond,[],2);
%     map_cond = roi_time(mode_idcs);
%     plot(map_ref,map_cond,...
%         'color',contrast_clrs(ii,:),...
%         'linewidth',1.5);
% end
%
% % plot identity line
% plot(xlim,ylim,'--k');