%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% generate positive control
if ~isfield(data,'FakeFR')
    toso2021_simulateControls;
end

%% neuron selection

% selected for being good examples of i2-modulation
if strcmpi(task_str,'duration')
    neurons2plot = fliplr([...
        393,215,72,526]);
elseif strcmpi(task_str,'intensity')
    neurons2plot = [...
        19,22,30,61,66,70,100,111,112,115,...
        166,238,243,260,344,408,410];
end
n_neurons2plot = numel(neurons2plot);

%% construct Si-aligned, Ti- & Ii-split psths
ti_padd = [-0,0];
sdf_gain = .75;
sdf_spacing = 1.05;

% figure initialization
fig = figure(figopt,...
    'position',[489,343,460,420],...
    'name',sprintf('fake_modulation_%s',contrast_str));

% axes initialization
sps = gobjects(2,1);
sps(1) = subplot(1,2,1);
sps(2) = subplot(1,2,2);
set(sps,...
    axesopt.default,...
    'plotboxaspectratiomode','auto',...
    'clipping','off',...
    'ylim',[0,n_neurons2plot+1]+[-1,1]*.05*(n_neurons2plot+1),...
    'ycolor','none');
xxtick = unique([ti_padd(1);-pre_s1_delay;0;t_set;t_set(end)+ti_padd(2)]);
xxtick(xxtick == 0) = abs(xxtick(xxtick == 0));
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,[0,1e3,t_set(t1_mode_idx)])) = {''};
set(sps(1),...
    'xlim',[0,max(t_set)]+ti_padd,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel);
xxtick = unique([-ti_padd(1);0;-t_set;-(t_set(end)+ti_padd(2))]);
xxtick(xxtick == 0) = abs(xxtick(xxtick == 0));
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,-[0,1e3,t_set(t1_mode_idx)])) = {''};
set(sps(2),...
    'xlim',sort(-([0,max(t_set)]+ti_padd)),...
    'xtick',xxtick,...
    'xticklabel',xxticklabel);
xlabel(sps(1),'Time since S_2 onset (ms)');
xlabel(sps(2),'Time since S_2 offset (ms)');
ylabel(sps(1),'Firing rate (Hz)');
ylabel(sps(2),'Firing rate (Hz)');

% iterate through neurons
for nn = 1 : n_neurons2plot
    progressreport(nn,n_neurons2plot,'parsing neural data');
    neuron_flags = data.NeuronNumb == neurons2plot(nn);
    
    % iterate through intensities
    for ii = n_i : -1 : 1
        i2_flags = i2 == i_set(ii);
        trial_flags = ...
            valid_flags & ...
            neuron_flags & ...
            i2_flags;
        n_flagged_trials = sum(trial_flags);
        if n_flagged_trials == 0
            continue;
        end

        % time selection
        time2plot = min(xlim(sps(1))) + psthbin : psthbin : max(xlim(sps(1)));

        % patch average activity
        s2on_mu = psths.s2on(:,neurons2plot(nn),ii)';
        nan_flags = isnan(s2on_mu);
        if ii == n_i
            s2_min = min(s2on_mu);
            s2_range = range(s2on_mu);
        end
        s2on_mu = (s2on_mu - s2_min) / s2_range - .5;
        s2on_mu = s2on_mu * sdf_gain + .5 + (nn - 1) * sdf_spacing;
        s2on_xpatch = time2plot(~nan_flags);
        s2on_ypatch = s2on_mu(~nan_flags);
        s2on_ypatch(end) = nan;
        patch_alpha = s2on_trials_throughtime(~nan_flags);
        patch_alpha = patch_alpha ./ max(patch_alpha) .* ...
            range(alphabounds_mu) + alphabounds_mu(1);
        alpha_levels = unique(patch_alpha,'stable');
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(1),...
                [s2on_xpatch(alpha_flags),nan],...
                [s2on_ypatch(alpha_flags),nan],0,...
                'edgealpha',alpha_levels(aa),...
                'edgecolor',i2_clrs(ii,:),...
                'facecolor','none',...
                'linewidth',1.5);
        end
        
        % time selection
        time2plot = min(xlim(sps(2))) + psthbin : psthbin : max(xlim(sps(2)));
        
        % patch average activity
        s2off_mu = psths.s2off(:,neurons2plot(nn),ii)';
        nan_flags = isnan(s2off_mu);
        s2off_mu = (s2off_mu - s2_min) / s2_range - .5;
        s2off_mu = s2off_mu * sdf_gain + .5 + (nn - 1) * sdf_spacing;
        s2off_xpatch = time2plot(~nan_flags);
        s2off_ypatch = s2off_mu(~nan_flags);
        s2off_ypatch(end) = nan;
        patch_alpha = s2off_trials_throughtime(~nan_flags);
        patch_alpha = patch_alpha ./ max(patch_alpha) .* ...
            range(alphabounds_mu) + alphabounds_mu(1);
        alpha_levels = unique(patch_alpha,'stable');
        n_alpha_levels = numel(alpha_levels);
        for aa = 1 : n_alpha_levels
            alpha_flags = patch_alpha == alpha_levels(aa);
            patch(sps(2),...
                [s2off_xpatch(alpha_flags),nan],...
                [s2off_ypatch(alpha_flags),nan],0,...
                'edgealpha',alpha_levels(aa),...
                'edgecolor',i2_clrs(ii,:),...
                'facecolor','none',...
                'linewidth',1.5);
        end
    end
end

% iterate through durations
for tt = 1 : n_t
    multiplier = 1 - .015 * (tt - 1);
    
    % plot spike integration window
    plot(sps(1),[0,t_set(tt)],[1,1]*max(ylim(sps(1)))*multiplier,'-k',...
        'linewidth',3);
    plot(sps(2),[-t_set(tt),0],[1,1]*max(ylim(sps(2)))*multiplier,'-k',...
        'linewidth',3);
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end