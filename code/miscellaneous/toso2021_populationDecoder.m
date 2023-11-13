%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% run settings
n_runs = 1;

%% session selection criteria
min_neuron_count = 10;

%% bootstrap settings
n_boots = 0; % 1e3;

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% eventually it'd be nice to drop neurons and keep tabs
% on performance, what's the kosher way of doing this?
% do a bunch of runs changing drop order?
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% SHOULD ITERATE THROUGH SESSIONS FIRST, SINCE THIS SHOUÇD BE DONE
% WITHIN SESSION AND THEN AVERAGED ACROSS SESSIONS (?). aRBITRATE A
% MINIMUM NUMBER OF NEURONS?

%% decoder settings
glm_wins = t_set(1);
n_glm = numel(glm_wins);

%% preallocation
spkcounts = struct();
performance = struct();

%% iterate through runs
for rr = 1 : n_runs
    if n_runs > 1
        toso2021_simulateSpikes;
        spike_data_field = 'FakeFR';
    else
        spike_data_field = 'FR';
    end
    
    %% iterate through spike integration windows
    for gg = 1 : n_glm
        glm_win = glm_wins(gg);
        glm_str = sprintf('t%i',glm_win);
        
        %% ROI settings
        
        % initialization
        glm_roi = struct();
        glm_roi_lbl = struct();
        
        % roi definitions
        glm_roi.preInit = [-glm_win,0];
        glm_roi.postInit = [0,glm_win];
        glm_roi.preS1Onset = [-glm_win,0];
        glm_roi.postS1Onset = [0,glm_win];
        glm_roi.preS1Offset = [-glm_win,0];
        glm_roi.postS1Offset = [0,glm_win];
        glm_roi.preS2Onset = [-glm_win,0];
        glm_roi.postS2Onset = [0,glm_win];
        glm_roi.preS2Offset = [-glm_win,0];
        glm_roi.postS2Offset = [0,glm_win];
        glm_roi.preGoCue = [-glm_win,0];
        glm_roi.postGoCue = [0,glm_win];
        glm_roi.aroundReaction = glm_roi.postGoCue + glm_win;
        glm_roi.aroundChoice = glm_roi.postGoCue + glm_win * 2;
        
        % roi labels
        glm_roi_lbl.preInit = 'Pre-initiation';
        glm_roi_lbl.postInit = 'Post-initiation';
        glm_roi_lbl.preS1Onset = 'Pre-S1 onset';
        glm_roi_lbl.postS1Onset = 'Post-S1 onset';
        glm_roi_lbl.preS1Offset = 'Pre-S1 offset';
        glm_roi_lbl.postS1Offset = 'Post-S1 offset';
        glm_roi_lbl.preS2Onset = 'Pre-S2 onset';
        glm_roi_lbl.postS2Onset = 'Post-S2 onset';
        glm_roi_lbl.preS2Offset = 'Pre-S2 offset';
        glm_roi_lbl.postS2Offset = 'Post-S2 offset';
        glm_roi_lbl.preGoCue = 'Pre-go cue';
        glm_roi_lbl.postGoCue = 'Post-go cue';
        glm_roi_lbl.aroundReaction = '~Reaction';
        glm_roi_lbl.aroundChoice = '~Choice';
        
        % adjust if simulating
        if n_runs > 1
            all_fields = fieldnames(glm_roi);
            fields2rm = ~ismember(all_fields,...
                {'postS1Onset','preS1Offset','postS2Onset','preS2Offset'});
            glm_roi = rmfield(glm_roi,all_fields(fields2rm));
            glm_roi_lbl = rmfield(glm_roi_lbl,all_fields(fields2rm));
        end
        
        % epoch parsing
        epochs = fieldnames(glm_roi);
        n_epochs = numel(epochs);
        
        %% iterate through sessions
        for ss = 1 : n_total_sessions
            session_flags = session_idcs == ss;
            
            %% neuron count-based selection
            if session_neuron_count(ss) <= min_neuron_count
                continue;
            end
            
            %% preallocation
            for ii = 1 : n_epochs
                spkcounts.(epochs{ii}) = ...
                    nan(n_glm,session_neuron_count(ss),session_trial_count(ss));
            end
            
            %% neuron selection
            session_neurons = unique(data.NeuronNumb(...
                valid_flags & ...
                session_flags & ...
                ismember(data.NeuronNumb,flagged_neurons)));
            
            %% construct response

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
                spike_counts = data.(spike_data_field)(trial_flags,:)';
%                 spike_counts = data.SDF(trial_flags,:)';
                n_trials = sum(trial_flags);
                
                % around approach spike rates
                if isfield(glm_roi,'aroundInitMov')
                    alignment_onset = ...
                        repmat(pre_init_padding,n_trials,1);
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.aroundInitMov(1) & ...
                        padded_time < alignment_onset + glm_roi.aroundInitMov(2);
                    chunk_flags = alignment_flags;
                    spkcounts_aroundInitMov = spike_counts;
                    spkcounts_aroundInitMov(~alignment_flags') = nan;
                    spkcounts_aroundInitMov = ...
                        reshape(spkcounts_aroundInitMov(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.aroundInitMov(gg,nn,:) = nansum(spkcounts_aroundInitMov,2);
                end
                
                % pre initiation spike rates
                if isfield(glm_roi,'preInit')
                    alignment_onset = ...
                        repmat(pre_init_padding,n_trials,1);
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.preInit(1) & ...
                        padded_time < alignment_onset + glm_roi.preInit(2);
                    chunk_flags = alignment_flags;
                    spkcounts_preInit = spike_counts;
                    spkcounts_preInit(~alignment_flags') = nan;
                    spkcounts_preInit = ...
                        reshape(spkcounts_preInit(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.preInit(gg,nn,:) = nansum(spkcounts_preInit,2);
                end
                
                % post initiation spike rates
                if isfield(glm_roi,'postInit')
                    alignment_onset = ...
                        repmat(pre_init_padding,n_trials,1);
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.postInit(1) & ...
                        padded_time < alignment_onset + glm_roi.postInit(2);
                    chunk_flags = alignment_flags;
                    spkcounts_postInit = spike_counts;
                    spkcounts_postInit(~alignment_flags') = nan;
                    spkcounts_postInit = ...
                        reshape(spkcounts_postInit(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.postInit(gg,nn,:) = nansum(spkcounts_postInit,2);
                end
                
                % pre-S1 onset spike rates
                if isfield(glm_roi,'preS1Onset')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags);
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.preS1Onset(1) & ...
                        padded_time < alignment_onset + glm_roi.preS1Onset(2);
                    chunk_flags = alignment_flags;
                    spkcounts_preS1Onset = spike_counts;
                    spkcounts_preS1Onset(~alignment_flags') = nan;
                    spkcounts_preS1Onset = ...
                        reshape(spkcounts_preS1Onset(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.preS1Onset(gg,nn,:) = nansum(spkcounts_preS1Onset,2);
                end
                
                % post-S1 onset spike rates
                alignment_onset = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags);
                alignment_flags = ...
                    padded_time >= alignment_onset + glm_roi.postS1Onset(1) & ...
                    padded_time < alignment_onset + t1(trial_flags);
                chunk_flags = ...
                    padded_time >= alignment_onset + glm_roi.postS1Onset(1) & ...
                    padded_time < alignment_onset + glm_roi.postS1Onset(2);
                spkcounts_postS1Onset = spike_counts;
                spkcounts_postS1Onset(~alignment_flags') = nan;
                spkcounts_postS1Onset = ...
                    reshape(spkcounts_postS1Onset(chunk_flags'),[glm_win,n_trials])';
                
                % pre-S1 offset spike rates
                alignment_onset = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags);
                alignment_flags = ...
                    padded_time >= alignment_onset - t1(trial_flags) & ...
                    padded_time < alignment_onset + glm_roi.preS1Offset(2);
                chunk_flags = ...
                    padded_time >= alignment_onset + glm_roi.preS1Offset(1) & ...
                    padded_time < alignment_onset + glm_roi.preS1Offset(2);
                spkcounts_preS1Offset = spike_counts;
                spkcounts_preS1Offset(~alignment_flags') = nan;
                spkcounts_preS1Offset = ...
                    reshape(spkcounts_preS1Offset(chunk_flags'),[glm_win,n_trials])';
                
                % post-S1 offset spike rates
                if isfield(glm_roi,'postS1Offset')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags) + ...
                        t1(trial_flags);
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.postS1Offset(1) & ...
                        padded_time < alignment_onset + glm_roi.postS1Offset(2);
                    chunk_flags = alignment_flags;
                    spkcounts_postS1Offset = spike_counts;
                    spkcounts_postS1Offset(~alignment_flags') = nan;
                    spkcounts_postS1Offset = ...
                        reshape(spkcounts_postS1Offset(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.postS1Offset(gg,nn,:) = nansum(spkcounts_postS1Offset,2);
                end
                
                % pre-S2 onset spike rates
                if isfield(glm_roi,'preS2Onset')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags) + ...
                        t1(trial_flags) + ...
                        isi;
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.preS2Onset(1) & ...
                        padded_time < alignment_onset + glm_roi.preS2Onset(2);
                    chunk_flags = alignment_flags;
                    spkcounts_preS2Onset = spike_counts;
                    spkcounts_preS2Onset(~alignment_flags') = nan;
                    spkcounts_preS2Onset = ...
                        reshape(spkcounts_preS2Onset(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.preS2Onset(gg,nn,:) = nansum(spkcounts_preS2Onset,2);
                end
                
                % post-S2 onset spike rates
                alignment_onset = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags) + ...
                    isi;
                alignment_flags = ...
                    padded_time >= alignment_onset + glm_roi.postS2Onset(1) & ...
                    padded_time < alignment_onset + t2(trial_flags);
                chunk_flags = ...
                    padded_time >= alignment_onset + glm_roi.postS2Onset(1) & ...
                    padded_time < alignment_onset + glm_roi.postS2Onset(2);
                spkcounts_postS2Onset = spike_counts;
                spkcounts_postS2Onset(~alignment_flags') = nan;
                spkcounts_postS2Onset = ...
                    reshape(spkcounts_postS2Onset(chunk_flags'),[glm_win,n_trials])';
                
                % pre-S2 offset spike rates
                alignment_onset = ...
                    pre_init_padding + ...
                    pre_s1_delay(trial_flags) + ...
                    t1(trial_flags) + ...
                    isi + ...
                    t2(trial_flags);
                alignment_flags = ...
                    padded_time >= alignment_onset - t2(trial_flags) & ...
                    padded_time < alignment_onset + glm_roi.preS2Offset(2);
                chunk_flags = ...
                    padded_time >= alignment_onset + glm_roi.preS2Offset(1) & ...
                    padded_time < alignment_onset + glm_roi.preS2Offset(2);
                spkcounts_preS2Offset = spike_counts;
                spkcounts_preS2Offset(~alignment_flags') = nan;
                spkcounts_preS2Offset = ...
                    reshape(spkcounts_preS2Offset(chunk_flags'),[glm_win,n_trials])';
                
                % post-S2 offset spike rates
                if isfield(glm_roi,'postS2Offset')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags) + ...
                        t1(trial_flags) + ...
                        isi + ...
                        t2(trial_flags);
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.postS2Offset(1) & ...
                        padded_time < alignment_onset + glm_roi.postS2Offset(2);
                    chunk_flags = alignment_flags;
                    spkcounts_postS2Offset = spike_counts;
                    spkcounts_postS2Offset(~alignment_flags') = nan;
                    spkcounts_postS2Offset = ...
                        reshape(spkcounts_postS2Offset(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.postS2Offset(gg,nn,:) = nansum(spkcounts_postS2Offset,2);
                end
                
                % pre-go spike rates
                if isfield(glm_roi,'preGoCue')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags) + ...
                        t1(trial_flags) + ...
                        isi + ...
                        t2(trial_flags) + ...
                        post_s2_delay;
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.preGoCue(1) & ...
                        padded_time < alignment_onset + glm_roi.preGoCue(2);
                    chunk_flags = alignment_flags;
                    spkcounts_preGoCue = spike_counts;
                    spkcounts_preGoCue(~alignment_flags') = nan;
                    spkcounts_preGoCue = ...
                        reshape(spkcounts_preGoCue(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.preGoCue(gg,nn,:) = nansum(spkcounts_preGoCue,2);
                end
                
                % post-go spike rates
                if isfield(glm_roi,'postGoCue')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags) + ...
                        t1(trial_flags) + ...
                        isi + ...
                        t2(trial_flags) + ...
                        post_s2_delay;
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.postGoCue(1) & ...
                        padded_time < alignment_onset + glm_roi.postGoCue(2);
                    chunk_flags = alignment_flags;
                    spkcounts_postGoCue = spike_counts;
                    spkcounts_postGoCue(~alignment_flags') = nan;
                    spkcounts_postGoCue = ...
                        reshape(spkcounts_postGoCue(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.postGoCue(gg,nn,:) = nansum(spkcounts_postGoCue,2);
                end
                
                % around choice spike rates
                if isfield(glm_roi,'aroundReaction')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags) + ...
                        t1(trial_flags) + ...
                        isi + ...
                        t2(trial_flags) + ...
                        post_s2_delay;
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.aroundReaction(1) & ...
                        padded_time < alignment_onset + glm_roi.aroundReaction(2);
                    chunk_flags = alignment_flags;
                    spkcounts_aroundReaction = spike_counts;
                    spkcounts_aroundReaction(~alignment_flags') = nan;
                    spkcounts_aroundReaction = ...
                        reshape(spkcounts_aroundReaction(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.aroundReaction(gg,nn,:) = nansum(spkcounts_aroundReaction,2);
                end
                
                % around choice spike rates
                if isfield(glm_roi,'aroundChoice')
                    alignment_onset = ...
                        pre_init_padding + ...
                        pre_s1_delay(trial_flags) + ...
                        t1(trial_flags) + ...
                        isi + ...
                        t2(trial_flags) + ...
                        post_s2_delay;
                    alignment_flags = ...
                        padded_time >= alignment_onset + glm_roi.aroundChoice(1) & ...
                        padded_time < alignment_onset + glm_roi.aroundChoice(2);
                    chunk_flags = alignment_flags;
                    spkcounts_aroundChoice = spike_counts;
                    spkcounts_aroundChoice(~alignment_flags') = nan;
                    spkcounts_aroundChoice = ...
                        reshape(spkcounts_aroundChoice(chunk_flags'),[glm_win,n_trials])';
                    spkcounts.aroundChoice(gg,nn,:) = nansum(spkcounts_aroundChoice,2);
                end
                
                % store average spike rates
                spkcounts.postS1Onset(gg,nn,:) = nansum(spkcounts_postS1Onset,2);
                spkcounts.preS1Offset(gg,nn,:) = nansum(spkcounts_preS1Offset,2);
                spkcounts.postS2Onset(gg,nn,:) = nansum(spkcounts_postS2Onset,2);
                spkcounts.preS2Offset(gg,nn,:) = nansum(spkcounts_preS2Offset,2);
            end
            
            %% LDA-based population decoder (1 vs. all)
            
            % trial selection
            trial_flags = ...
                valid_flags & ...
                unique_flags & ...
                session_flags;
            
            % response selection
            response_table = table(t1,t2,i1,i2);
            response_table = response_table(trial_flags,:);
            response_names = response_table.Properties.VariableNames;
            n_responses = size(response_table,2);

            % iterate through responses
            for ii = 1 : n_responses
                response_str = response_names{ii};
                response = response_table.(response_str);
                response_set = eval([response_str,'_set']);
                response_clrs = eval([response_str,'_clrs']);
                n_response_values = numel(response_set);
                
                % iterate through epochs
                for ee = 1 : n_epochs
                    epoch = epochs{ee};
                    
                    % design matrix
                    X = squeeze(spkcounts.(epoch)(gg,:,:))';
                    
%                     if strcmpi(epoch,'postS2Onset');
%                         figure;
%                         [~,idcs] = sort(response);
%                         imagesc([X(idcs,:),response(idcs)/10]);
%                         a=1
%                     end
                    
%                     figure('windowstyle','docked');
%                     imagesc([X,response]);
%                     a=1

%                     figure; hold on;
%                     title(epoch);
                    
                    % iterate through response values
                    for jj = 1 % : n_response_values
%                         plot(nanmean(X(response == response_set(jj),:)),...
%                             'color',response_clrs(jj,:),...
%                             'marker','.',...
%                             'markersize',25);
%                         continue;
                        
                        % response variable
                        y = response == response_set(jj);
                        
                        % iterate through trials
                        for tt = 1 : session_trial_count(ss)
                            progressreport(tt,session_trial_count(ss),...
                                sprintf('cross-validating (%s, %s = %i)',...
                                epoch,response_str,response_set(jj)));
                            
                            % handle leave-one-out cross-validation
                            train_flags = (1 : session_trial_count(ss)) ~= tt;
                            train_idcs = find(train_flags);
                            shuffled_idcs = ...
                                train_idcs(randperm(session_trial_count(ss)-1));
                            
                            try
                            % linear discriminant analysis
%                             mdl = fitcdiscr(X(train_flags,:),y(train_flags),...
%                                 'discrimtype','linear');
                            mdl = fitcecoc(...
                                X(train_flags,:),response(train_flags),...
                                'learners','discriminant',...
                                'coding','onevsall');
                            
                            % prediction with test trial
%                             performance.(epoch)(gg,jj,tt) = ...
%                                 mdl.predict(X(tt,:)) == y(tt);
                            performance.(epoch)(gg,jj,tt) = ...
                                mdl.predict(X(tt,:)) == response(tt);
                            catch
                            end
                        end
                    end
                end
                
                %% plot decoding performance
                
                % bar width settings
                epoch_span = 1 / 3;
                barwidth = epoch_span / n_responses;
                
                % figure initialization
                fig = figure(figopt,...
                    'name',sprintf('population_decoder_%s_s%i',response_str,ss),...
                    'windowstyle','docked',...
                    'color',[1,1,1]*1);
                
                % axes initialization
                axes(axesopt.default,...
                    'plotboxaspectratio',[1,1,1],...
                    'color','none',...
                    'xlim',[1,n_epochs]+[-1,1]*.05*n_epochs,...
                    'xtick',1:n_epochs,...
                    'xticklabel',num2cell(1:n_epochs),...
                    'xticklabelrotation',90,...
                    'ylim',[0,1]+[-1,1]*.05,...
                    'ytick',linspace(0,1,5),...
                    'yticklabel',num2cell(linspace(0,1,5)*100),...
                    'xcolor','k',...
                    'clipping','off',...
                    'layer','bottom');
                xlabel('Task epoch');
                ylabel('Decoder performance (%)');
                
                % reference lines
                plot(xlim,[1,1]*0,'-k',...
                    'linewidth',1.5);
                
                % pseudo-legend (stimulus epochs)
                s1epoch_idcs = ...
                    [find(ismember(epochs,'postS1Onset')),...
                    find(ismember(epochs,'preS1Offset'))];
                plot(s1epoch_idcs+[-1,1]*.5,[1,1],...
                    'linestyle','-',...
                    'linewidth',3,...
                    'color',stim_clrs(1,:));
                text(mean(s1epoch_idcs),1.025,'S1',...
                    'color',stim_clrs(1,:),...
                    'fontsize',12,...
                    'horizontalalignment','center',...
                    'verticalalignment','bottom');
                s2epoch_idcs = ...
                    [find(ismember(epochs,'postS2Onset')),...
                    find(ismember(epochs,'preS2Offset'))];
                plot(s2epoch_idcs+[-1,1]*.5,[1,1],...
                    'linestyle','-',...
                    'linewidth',3,...
                    'color',stim_clrs(2,:));
                text(mean(s2epoch_idcs),1.025,'S2',...
                    'color',stim_clrs(2,:),...
                    'fontsize',12,...
                    'horizontalalignment','center',...
                    'verticalalignment','bottom');

                % iterate through epochs
                for ee = 1 : n_epochs
                    epoch = epochs{ee};
                    
                    % epoch delimeter
                    p = plot([1,1]*ee+1/2,[0,1],...
                        'color','k',...
                        'linestyle','--');
                    uistack(p,'bottom');
                    
                    % pseudo x-tick label
                    text(ee,-.1,glm_roi_lbl.(epoch),...
                        'horizontalalignment','right',...
                        'rotation',90,...
                        'color','k',...
                        'fontsize',10);
                    
                    % compute horizontal offsets
                    x_offsets = ee + (0 : n_epochs - 1) * barwidth + ...
                        barwidth / 2 - epoch_span / 2;
                    x = x_offsets(ee);
                    
                    % iterate through response values
                    for jj = 1 : n_response_values
                        try
                        plot(ee,nanmean(performance.(epoch)(gg,jj,:)),...
                            'marker','o',...
                            'markersize',7.5,...
                            'markerfacecolor',response_clrs(jj,:),...
                            'markeredgecolor','k',...
                            'linewidth',1.5);
                        catch 
                        end
                    end
                    continue;
                end
                
                % save figure
                if want2save
                    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
                    print(fig,svg_file,'-dsvg','-painters');
                end
            end
        end
    end
end