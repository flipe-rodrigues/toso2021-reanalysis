%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

for ss = 1 : n_total_sessions
    session_flags = session_idcs == ss;
    
    %% construct psychophysical triple
    
    % normalize stimulus range
    norm_stimuli = (stimuli - min(stimuli)) / range(stimuli);
    normstim_set = unique(norm_stimuli(valid_flags));
    
    % preallocation
    psycurves = struct();
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        contrast_flags = contrasts == contrast_set(kk);
        
        % iterate through stimuli
        for ii = 1 : n_stimuli
            stimulus_flags = stimuli == stim_set(ii);
            trial_flags = ...
                session_flags & ...
                valid_flags & ...
                unique_flags & ...
                contrast_flags & ...
                stimulus_flags;
            psycurves(kk).x(ii,1) = normstim_set(ii);
            psycurves(kk).y(ii,1) = sum(choice(trial_flags));
            psycurves(kk).n(ii,1) = sum(trial_flags);
            psycurves(kk).err(ii,1) = ...
                std(choice(trial_flags)) / sqrt(sum(trial_flags));
        end
    end
    
    %% fit psychometric function
    
    % psychometric fit settings
    psyopt.fit = struct();
    psyopt.fit.expType = 'YesNo';
    psyopt.fit.sigmoidName = 'gauss';
    psyopt.fit.estimateType = 'MAP';
    psyopt.fit.confP = [.95,.9,.68];
    psyopt.fit.borders = [0,1; 0,1; 0,.25; 0,.25; 0,0];
    psyopt.fit.fixedPars = [nan,nan,nan,nan,0];
    psyopt.fit.stepN = [100,100,40,40,20];
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        
        % fit psychometric curve
        psycurves(kk).fit = ...
            psignifit([psycurves(kk).x,psycurves(kk).y,psycurves(kk).n],psyopt.fit);
    end
    
    % fit big psychometric curve
    % bigpsy.fit = ...
    %     psignifit([bigpsy.x,bigpsy.y,bigpsy.n],psyopt.fit);
    
    %% plot phychometric function
    
    % stimulus-specific axes properties
    axesopt.stimulus.xlim = ...
        ([normstim_set(1),normstim_set(end)] +  [-1,1] * .05);
    axesopt.stimulus.xtick = normstim_set;
    axesopt.stimulus.xticklabel = num2cell(round(stim_set,2));
    ticks2delete = ...
        ~ismember(axesopt.stimulus.xtick,...
        [min(axesopt.stimulus.xtick),max(axesopt.stimulus.xtick)]);
    axesopt.stimulus.xticklabel(ticks2delete) = {''};
    
    % psychometric plot settings
    psyopt.plot = struct();
    psyopt.plot.linewidth = 1.5;
    psyopt.plot.marker = 'o';
    psyopt.plot.markersize = 8.5;
    psyopt.plot.plotdata = true;
    psyopt.plot.gradeclrs = false;
    psyopt.plot.patchci = false;
    psyopt.plot.normalizemarkersize = true;
    
    % figure initialization
    fig = figure(figopt,...
        'windowstyle','docked',...
        'name',sprintf('psychometric_curves_%s',contrast_str));
    
    % axes initialization
    axes(...
        axesopt.default,...
        axesopt.stimulus,...
        axesopt.psycurve);
    title(sprintf('animal %i, session %i, %i neurons',...
        unique(subjects(session_flags&valid_flags&unique_flags)),ss,...
        numel(unique(flagged_neurons(session_flags&valid_flags&unique_flags)))));
    xlabel('T_2 (ms)');
    ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));
    
    % graphical object preallocation
    p = gobjects(n_contrasts,1);
    
    % reference lines
    plot([1,1]*median(normstim_set),ylim,':k');
    plot(xlim,[1,1]*.5,':k');
    
    % iterate through contrasts
    for kk = 1 : n_contrasts
        
        % plot psychometric curve
        psyopt.plot.datafaceclr = contrast_clrs(kk,:);
        psyopt.plot.overallvisibility = 'off';
        psyopt.plot.plotfit = sum(psycurves(kk).n ~= 0) >= 2;
        p(kk) = plotpsy(psycurves(kk),psycurves(kk).fit,psyopt.plot);
    end
    
    % legend
    if iscategorical(contrasts)
        leg_str = categories(contrast_set);
    else
        leg_str = cellfun(@(x,y)sprintf('%s = %i %s',x,y,contrast_units),...
            repmat({contrast_lbl},n_contrasts,1),num2cell(contrast_set),...
            'uniformoutput',false);
    end
    legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
        'edgecolor','k',...
        'position',[0.085,0.66,.27,.2],...
        'box','on');
end