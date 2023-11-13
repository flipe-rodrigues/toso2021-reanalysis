%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% session selection criteria
min_neuron_count = 9;

%% neurometric curve settings
spk_integration_win = t_set(1);
n_integration_bins = spk_integration_win / psthbin;

%% preallocation
neurocurves = struct();

%% neurometric curves
session_counter = 0;

% iterate through sessions
for ss = 1 : n_total_sessions
    session_flags = session_idcs == ss;
    
    % neuron count-based selection
    if session_neuron_count(ss) < min_neuron_count || ...
            session_subject(ss) == 0
        continue;
    end
    
    %% construct population vector
    
    % preallocation
    population_vector = nan(...
        session_neuron_count(ss),...
        session_trial_count(ss));
    
    % neuron selection
    session_neurons = unique(data.NeuronNumb(...
        valid_flags & ...
        session_flags & ...
        ismember(data.NeuronNumb,flagged_neurons)));
    
    % increment session index
    session_counter = session_counter + 1;
    
    % iterate through neurons
    for nn = 1 : session_neuron_count(ss)
        progressreport(nn,session_neuron_count(ss),...
            'fetching spike counts');
        neuron_flags = data.NeuronNumb == session_neurons(nn);
        
        % flag trials for the current condition
        trial_flags = ...
            valid_flags & ...
            session_flags & ...
            neuron_flags;
        if sum(trial_flags) == 0
            continue;
        end
        
        % fetch spike counts & compute spike rates
        spike_counts = data.SDF(trial_flags,:)';
        n_trials = sum(trial_flags);
        
        % S2-offset-aligned spike counts
        alignment_onset = ...
            pre_init_padding + ...
            pre_s1_delay(trial_flags) + ...
            t1(trial_flags) + ...
            isi + ...
            t2(trial_flags);
        alignment_flags = ...
            padded_time >= alignment_onset - spk_integration_win & ...
            padded_time < alignment_onset;
        chunk_flags = alignment_flags;
        aligned_spkcounts = spike_counts;
        aligned_spkcounts(~alignment_flags') = nan;
        aligned_spkcounts = reshape(aligned_spkcounts(chunk_flags'),...
            [n_integration_bins,n_trials])';
        
        % store spike count
        population_vector(nn,:) = nanmean(aligned_spkcounts,2);
    end
    
    %% LDA-based population decoder
    
    % trial selection
    trial_flags = ...
        valid_flags & ...
        unique_flags & ...
        session_flags;
    
    % preallocation
    neural_choices = nan(session_trial_count(ss),1);
    
    % design matrix
    X = population_vector';
%     X = (X - nanmean(X,1)) ./ std(X,0,1);
    
    % response variable
%     y = choice(trial_flags);
    
%     threshold = median(s1(valid_flags));
%     ambiguous_flags = stimuli(trial_flags) == threshold;
%     y = stimuli(trial_flags) > threshold;
%     y(ambiguous_flags) = rand(sum(ambiguous_flags),1) > .5;
    
%     weighted_stimuli = s2(trial_flags) * w2 + s1(trial_flags) * w1;
%     y = weighted_stimuli > 200;
%     
    y = s2(trial_flags) > s1(trial_flags);

    % iterate through trials
    for tt = 1 : session_trial_count(ss)
        progressreport(tt,session_trial_count(ss),...
            sprintf('cross-validating (session %i)',ss));
        
        % handle leave-one-out cross-validation
        xval_flags = (1 : session_trial_count(ss))' ~= tt;
        train_flags = ...
            ...t1(trial_flags) == t_set(t1_mode_idx) & ...
            ...i1(trial_flags) == i_set(i1_mode_idx) & ...
            contrasts(trial_flags) == contrast_set(contrast_mode_idx) & ...
            xval_flags;
        
        % response variable
%         y(ambiguous_flags) = rand(sum(ambiguous_flags),1) > .5;
        
        % flag invalid neurons
        within_class_var = cell2mat(...
            arrayfun(@(x)nanvar(X(train_flags & y==x,:)),unique(y),...
            'uniformoutput',false));
        invalid_flags = any(...
            within_class_var == 0 | ...
            isnan(within_class_var));
        
        % linear discriminant analysis
        mdl = fitcdiscr(...
            X(train_flags,~invalid_flags),y(train_flags),...
            'discrimtype','linear');
        %         mdl = fitcsvm(X(train_flags,~invalid_flags),y(train_flags));
        
        % neural "judgments"
        neural_choices(tt) = mdl.predict(X(tt,~invalid_flags));
    end
    
    %% construct neurophysical triples
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        contrast_flags = contrasts(trial_flags) == contrast_set(kk);
        
        % iterate through stimuli
        for ii = 1 : n_stimuli
            stimulus_flags = stimuli(trial_flags) == stim_set(ii);
            flags = ...
                contrast_flags & ...
                stimulus_flags;
            
            % construct neurophysical triple
            neurocurves(kk).x(ii,session_counter) = normstim_set(ii);
            neurocurves(kk).y(ii,session_counter) = ...
                nansum(neural_choices(flags),[1,2]);
            neurocurves(kk).n(ii,session_counter) = nansum(flags);
            neurocurves(kk).err(ii,session_counter) = ...
                std(neural_choices(flags),0,[1,2]) / sqrt(sum(flags));
        end
    end
end

%% fit neurometric functions

% preallocation
neurocurve_pools = struct();

% neurometric fit settings
psyopt.fit = struct();
psyopt.fit.expType = 'YesNo';
psyopt.fit.sigmoidName = 'gauss';
psyopt.fit.estimateType = 'MAP';
psyopt.fit.confP = [.95,.9,.68];
psyopt.fit.borders = [0,1; 0,1; 0,.4; 0,.4; 0,0];
psyopt.fit.fixedPars = [nan,nan,nan,nan,0];
psyopt.fit.stepN = [100,100,40,40,20];

% iterate through contrasts
for kk = 1 : n_contrasts
    
    % pool neural choices across sessions
    neurocurve_pools(kk).x = unique(neurocurves(kk).x);
    neurocurve_pools(kk).p = ...
        nanmean(neurocurves(kk).y ./ neurocurves(kk).n,2);
    neurocurve_pools(kk).y = floor(...
        neurocurve_pools(kk).p .* nansum(neurocurves(kk).n,2));
    neurocurve_pools(kk).n = nansum(neurocurves(kk).n,2);
    neurocurve_pools(kk).err = [1,1] .* ...
        nanstd(neurocurves(kk).y ./ neurocurves(kk).n,0,2) / sqrt(session_counter);
    
    % neurophysical triple
    neuro_triple = [...
        neurocurve_pools(kk).x,...
        neurocurve_pools(kk).y,...
        neurocurve_pools(kk).n...
        ];
    
    % fit neurometric curve
    neurocurve_pools(kk).fit = psignifit(neuro_triple,psyopt.fit);
end

% % average & error function selection
% avgfun = @(x,d) nansum(x,d);
% errfun = @(x,d) std(x,0,d) / sqrt(size(x,2)) .* [1,1];
%
% % iterate through contrasts
% for kk = 1 : n_contrasts
%
%     % pool neural choices across sessions
%     neurocurve_pools(kk).x = unique(neurocurves(kk).x);
%     neurocurve_pools(kk).y = avgfun(neurocurves(kk).y,2);
%     neurocurve_pools(kk).n = avgfun(neurocurves(kk).n,2);
%     neurocurve_pools(kk).err = ....
%         errfun(neurocurves(kk).y ./ neurocurves(kk).n,2);
%
%     % neurophysical triple
%     neuro_triple = [...
%         neurocurve_pools(kk).x,...
%         neurocurve_pools(kk).y,...
%         neurocurve_pools(kk).n...
%         ];
%
%     % fit neurometric curve
%     neurocurve_pools(kk).fit = psignifit(neuro_triple,psyopt.fit);
% end

%% plot neurometric functions

% figure initialization
fig = figure(figopt,...
    'name',sprintf('neurometric_curves_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);
xlabel('T_2 (ms)');
ylabel('P(T_2 > T_1)');

% neurometric plot settings
psyopt.plot = struct();
psyopt.plot.linewidth = 1.5;
psyopt.plot.marker = 's';
psyopt.plot.markersize = 9.5;
psyopt.plot.plotdata = true;
psyopt.plot.gradeclrs = false;
psyopt.plot.patchci = false;
psyopt.plot.overallvisibility = 'off';
psyopt.plot.normalizemarkersize = true;
psyopt.plot.plotfit = true;

% graphical object preallocation
p = gobjects(n_contrasts,1);

% reference lines
plot([1,1]*median(normstim_set),ylim,':k');
plot(xlim,[1,1]*.5,':k');

% iterate through contrasts
for kk = 1 : n_contrasts
    
    % plot neurometric curve
    psyopt.plot.datafaceclr = contrast_clrs(kk,:);
    psyopt.plot.overallvisibility = 'off';
    psyopt.plot.plotfit = sum(neurocurve_pools(kk).n ~= 0) >= 2;
    p(kk) = plotpsy(neurocurve_pools(kk),neurocurve_pools(kk).fit,psyopt.plot);
end

%% inset with delta P(T2 > T1)
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.inset.se,...
    'ylim',[-.25,.25]+[-1,1]*.05*.75,...
    'ytick',[-.25,0,.25],...
    'clipping','off');
ylabel('\DeltaP(T_2 > T_1)',...
    'rotation',-90,...
    'verticalalignment','bottom');

% categorical boundary
plot([1,1]*median(normstim_set),ylim,...
    'color','k',...
    'linestyle',':');

% zero level
plot(xlim,[0,0],...
    'color','k',...
    'linestyle',':');

% iterate through contrasts
for kk = 1 : n_contrasts
    if kk == contrast_mode_idx
        continue;
    end
    
    % plot delta P(T2 > T1)
    p_ctrl = neurocurve_pools(contrast_mode_idx).y ./ ...
        neurocurve_pools(contrast_mode_idx).n;
    p_cond = neurocurve_pools(kk).y ./ neurocurve_pools(kk).n;
    delta_p = p_cond - p_ctrl;
    err = sqrt((neurocurve_pools(kk).err .^ 2) + (neurocurve_pools(kk).err .^ 2));
    errorbar(normstim_set,delta_p,err(:,1),err(:,2),...
        'color','k',...
        'marker',psyopt.plot.marker,...
        'markersize',psyopt.plot.markersize-2.5,...
        'markerfacecolor',contrast_clrs(kk,:),...
        'linewidth',psyopt.plot.linewidth,...
        'capsize',0);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end
return;

%% fit & plot individual sessions
for ss = 1 : session_counter
    progressreport(ss,session_counter,'individual sessions');
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        
        ss_neurocurve(kk).x = neurocurves(kk).x(:,ss);
        ss_neurocurve(kk).y = neurocurves(kk).y(:,ss);
        ss_neurocurve(kk).n = neurocurves(kk).n(:,ss);
            ss_neurocurve(kk).err = ...
                std(ss_neurocurve(kk).y./ss_neurocurve(kk).n) ./ ...
                sqrt(ss_neurocurve(kk).n);
        
        % neurophysical triple
        neuro_triple = [...
            ss_neurocurve(kk).x,...
            ss_neurocurve(kk).y,...
            ss_neurocurve(kk).n...
            ];
        
        % fit neurometric curve
        ss_neurocurve(kk).fit = psignifit(neuro_triple,psyopt.fit);
    end
    
    % figure initialization
    fig = figure(figopt,...
        'windowstyle','docked',...
        'name',sprintf('neurometric_curves_%s_%i',contrast_str,ss));
    
    % axes initialization
    axes(...
        axesopt.default,...
        axesopt.stimulus,...
        axesopt.psycurve);
    xlabel('T_2 (ms)');
    ylabel('P(T_2 > T_1)');
    
    % reference lines
    plot([1,1]*median(normstim_set),ylim,':k');
    plot(xlim,[1,1]*.5,':k');
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        
        % plot neurometric curve
        psyopt.plot.datafaceclr = contrast_clrs(kk,:);
        psyopt.plot.overallvisibility = 'off';
        psyopt.plot.plotfit = true;
        psyopt.plot.normalizemarkersize = false;
        plotpsy(ss_neurocurve(kk),ss_neurocurve(kk).fit,psyopt.plot);
    end
end