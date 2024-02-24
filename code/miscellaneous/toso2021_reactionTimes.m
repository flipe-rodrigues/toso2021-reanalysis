%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% need to parse things first...
rt_i1 = data.Rts(:,3);
rt_t1 = data.Rts(:,4);
rt_inter_t1t2_edlay = data.Rts(:,5);
rt_i2 = data.Rts(:,6);
rt_t2 = data.Rts(:,7);
rt_ch = data.Rts(:,8);
rt = data.Rts(:,9);

%%
rt_i1(rt_i1 == 79) = 80;
rt_i1(rt_i1 == 149) = 147;
rt_i2(rt_i2 == 79) = 80;
rt_i2(rt_i2 == 149) = 147;
rt_valid_flags = ...
    ...data.Rts(:,1) ~= 1 & ...
    ...data.Rts(:,1) ~= 5 & ...
    rt_i1 ~= 0 & ...
    rt_i2 ~= 0 & ...
    rt_i2 ~= 53 & ...
    rt_i2 ~= 65 & ...
    rt_i2 ~= 120;

%% convert intensity from standard deviation to mean speed units
rt_i1 = round(rt_i1 .* sqrt(2 / pi));
rt_i2 = round(rt_i2 .* sqrt(2 / pi));

%% multiplicity
% rt_stimuli = round(rt_t2 * w2 + rt_t1 * w1);
% t2_minus_t1 = rt_stimuli(rt_valid_flags) > ...
%     nanmedian(unique((rt_stimuli(rt_valid_flags))));
% t2_plus_t1 = rt_t2(rt_valid_flags) + rt_t1(rt_valid_flags) > ...
%     nanmedian(unique((rt_t2(rt_valid_flags) + rt_t1(rt_valid_flags))));
% t2_alone = rt_t2(rt_valid_flags) > nanmedian(unique((rt_t2(rt_valid_flags))));
%
% ground_truth = rt_t2(rt_valid_flags) - rt_t1(rt_valid_flags) > 0;
%
% figure; hold on;
% win = 50;
% kernel = expkernel('mus',win,'binwidth',1);
% % plot(smooth(ground_truth == rt_ch(rt_valid_flags),win));
% % plot(smooth(t2_minus_t1 == rt_ch(rt_valid_flags),win));
% % plot(smooth(t2_plus_t1 == rt_ch(rt_valid_flags),win));
% % plot(smooth(t2_alone == rt_ch(rt_valid_flags),win));
% plot(conv(rt_ch(rt_valid_flags) == ground_truth,kernel.pdf,'valid'),...
%     'linewidth',1.5);
% plot(conv(t2_minus_t1 == ground_truth,kernel.pdf,'valid'));
% plot(conv(t2_plus_t1 == ground_truth,kernel.pdf,'valid'));
% plot(conv(t2_alone == ground_truth,kernel.pdf,'valid'));
%
% ylabel('Proportion')
% xlabel('Trial #')
% xlim([win,5000])
% ylim([0,1])
%
% legend({...
%     'choices == T_2 - T_1',...
%     'w_2 \times T_2 + w_1 \times T_1 == T_2 - T_1',...
%     'T_2 + T_1 == T_2 - T_1',...
%     'T_2 > med(T_2) == T_2 - T_1'})

%%
% contrast settings
rt_contrast = rt_i2;
rt_contrast_set = unique(rt_contrast(rt_valid_flags));
rt_n_contrasts = numel(rt_contrast_set);
rt_contrast_clrs = copper(rt_n_contrasts);
rt_stim = rt_t2; % - rt_t1;
rt_stim_set = unique(rt_stim(rt_valid_flags));
rt_n_stimuli = numel(rt_stim_set);

% normalize stimulus range
rt_norm_stimuli = (rt_stim - min(rt_stim)) / range(rt_stim);
rt_normstim_set = unique(rt_norm_stimuli(rt_valid_flags));

%% compute average reaction times

% preallocation
rts = struct();

% iterate through contrasts
for kk = 1 : rt_n_contrasts
    rt_contrast_flags = rt_contrast == rt_contrast_set(kk);
    
    % iterate through stimuli
    for ii = 1 : rt_n_stimuli
        rt_stim_flags = rt_stim == rt_stim_set(ii);
        
        % iterate through correctness
        for cc = [0,1]
            rt_choice_flags = rt_ch == cc; % (rt_ch == (rt_stim > nanmedian(rt_stim))) == cc;
            trial_flags = ...
                rt_valid_flags & ...
                rt_contrast_flags & ...
                rt_stim_flags & ...
                rt_choice_flags;
            rts(kk,cc+1).x(ii,1) = rt_normstim_set(ii);
            rts(kk,cc+1).y(ii,1) = nanmedian(rt(trial_flags));
            rts(kk,cc+1).n(ii,1) = sum(trial_flags);
            rts(kk,cc+1).err(ii,:) = ...
                ...[-1,1] * nanstd(rt(trial_flags)) / sqrt(sum(trial_flags));
                quantile(rt(trial_flags),[.25,.75]) - nanmedian(rt(trial_flags));
        end
    end
end

%% plot T2-split reaction times

% figure initialization
fig = figure(figopt,...
    'name',sprintf('reaction_times_t2split_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus);
title('Reaction times');
xlabel('T_2 (ms)');
ylabel('Reaction time (ms)');

% graphical object preallocation
p = gobjects(rt_n_contrasts,1);

% correctness markers
markers = {'s','o'};

% iterate through contrasts
for kk = 1 : rt_n_contrasts
    x_offset = (kk - ceil(rt_n_contrasts / 2)) * range(normstim_set) * .01;
    
    % iterate through correctness
    for cc = [0,1]
        
        % plot error bars
        errorbar(rts(kk,cc+1).x+x_offset,rts(kk,cc+1).y,...
            rts(kk,cc+1).err(:,1),rts(kk,cc+1).err(:,2),...
            'color',rt_contrast_clrs(kk,:),...
            'marker','none',...
            'linestyle','-',...
            'linewidth',1,...
            'capsize',2);
        
        % plot averages
        p(kk) = plot(rts(kk,cc+1).x+x_offset,rts(kk,cc+1).y,...
            'color',rt_contrast_clrs(kk,:),...
            'marker',markers{cc+1},...
            'markerfacecolor',rt_contrast_clrs(kk,:),...
            'markeredgecolor','k',...
            'linestyle','-',...
            'markersize',8.5,...
            'linewidth',1.5);
    end
end

% update axes limits
axis tight;
xlim(axesopt.stimulus.xlim);
ylim(ylim+[-1,1]*.05*range(ylim));

% reference lines
plot([1,1]*median(rt_normstim_set),ylim,':k');

% legend
leg_str = cellfun(@(x,y)sprintf('%s = %i %s',x,y,contrast_units),...
    repmat({contrast_lbl},rt_n_contrasts,1),num2cell(rt_contrast_set),...
    'uniformoutput',false);
legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
    'position',[0.085,0.66,.27,.2],...
    'box','on');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% plot reaction times

% figure initialization
fig = figure(figopt,...
    'name',sprintf('reaction_times_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    'xtick',unique(rt_i2(rt_valid_flags)));
title('Reaction times');
xlabel('I_2 (mm.s^{-1})');
ylabel('Reaction time (ms)');

% graphical object preallocation
p = gobjects(rt_n_contrasts,1);

% correctness markers
markers = {'s','o'};

% iterate through contrasts
for kk = 1 : rt_n_contrasts
    rt_contrast_flags = rt_contrast == rt_contrast_set(kk);
    
    % iterate through choices
    for cc = [0,1]
        rt_choice_flags = rt_ch == cc;
        trial_flags = ...
            rt_valid_flags & ...
            rt_contrast_flags & ...
            rt_choice_flags;
        
        rt_med = nanmedian(rt(trial_flags));
        rt_iqr = quantile(rt(trial_flags),[.25,.75]) - rt_med;
        
        x_offset = cc * 1.5;
        
        % plot error bars
        errorbar(rt_contrast_set(kk)+x_offset,rt_med,rt_iqr(1),rt_iqr(2),...
            'color',rt_contrast_clrs(kk,:),...
            'marker','none',...
            'linestyle','-',...
            'linewidth',1,...
            'capsize',2);
        
        % plot averages
        p(kk) = plot(rt_contrast_set(kk)+x_offset,rt_med,...
            'color',rt_contrast_clrs(kk,:),...
            'marker',markers{cc+1},...
            'markerfacecolor',rt_contrast_clrs(kk,:),...
            'markeredgecolor','k',...
            'linestyle','-',...
            'markersize',8.5+.25*(strcmpi(markers{cc+1},'s')),...
            'linewidth',1.5);
    end
end

% update axes limits
axis tight;
xlim(xlim+[-1,1]*.05*range(xlim));
ylim(ylim+[-1,1]*.05*range(ylim));

% legend
leg_str = cellfun(@(x,y)sprintf('%s = %i %s',x,y,contrast_units),...
    repmat({contrast_lbl},rt_n_contrasts,1),num2cell(rt_contrast_set),...
    'uniformoutput',false);
legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
    'position',[0.085,0.66,.27,.2],...
    'box','on');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end