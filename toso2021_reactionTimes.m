%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% need to parse things first...
rt_i1 = data.Rts(:,3);
rt_t1 = data.Rts(:,4);
rt_inter_t1t2_edlay = data.Rts(:,5);
rt_i2 = data.Rts(:,6);
rt_t2 = data.Rts(:,7);
rt_ch = data.Rts(:,8);
rt = data.Rts(:,9);

%%
rt_i2(rt_i2 == 79) = 80;
rt_i2(rt_i2 == 149) = 147;
rt_valid_flags = ...
    data.Rts(:,1) ~= 1 & ...
    data.Rts(:,1) ~= 3 & ...
    rt_i1 ~= 0 & ...
    rt_i2 ~= 0 & ...
    rt_i2 ~= 53 & ...
    rt_i2 ~= 65 & ...
    rt_i2 ~= 120;

%% multiplicity
rt_stimuli = round(rt_t2 * w2 + rt_t1 * w1);
t2_minus_t1 = rt_stimuli(rt_valid_flags) > ...
    nanmedian(unique((rt_stimuli(rt_valid_flags))));
t2_plus_t1 = rt_t2(rt_valid_flags) + rt_t1(rt_valid_flags) > ...
    nanmedian(unique((rt_t2(rt_valid_flags) + rt_t1(rt_valid_flags))));
t2_alone = rt_t2(rt_valid_flags) > nanmedian(unique((rt_t2(rt_valid_flags))));

ground_truth = rt_t2(rt_valid_flags) - rt_t1(rt_valid_flags) > 0;

figure; hold on;
win = 50;
kernel = expkernel('mus',win,'binwidth',1);
% plot(smooth(ground_truth == rt_ch(rt_valid_flags),win));
% plot(smooth(t2_minus_t1 == rt_ch(rt_valid_flags),win));
% plot(smooth(t2_plus_t1 == rt_ch(rt_valid_flags),win));
% plot(smooth(t2_alone == rt_ch(rt_valid_flags),win));
plot(conv(rt_ch(rt_valid_flags) == ground_truth,kernel.pdf,'valid'),...
    'linewidth',1.5);
plot(conv(t2_minus_t1 == ground_truth,kernel.pdf,'valid'));
plot(conv(t2_plus_t1 == ground_truth,kernel.pdf,'valid'));
plot(conv(t2_alone == ground_truth,kernel.pdf,'valid'));

ylabel('Proportion')
xlabel('Trial #')
xlim([win,5000])
ylim([0,1])

legend({...
    'choices == T_2 - T_1',...
    'w_2 \times T_2 + w_1 \times T_1 == T_2 - T_1',...
    'T_2 + T_1 == T_2 - T_1',...
    'T_2 > med(T_2) == T_2 - T_1'})

%%
rt_contrast = eval(['rt_',contrast_str]);
rt_stim = rt_t2 - rt_t1;

%% compute average reaction times

% preallocation
rts = struct();

% iterate through contrasts
for kk = 1 : n_contrasts
    contrast_flags = rt_contrast == contrast_set(kk);
    
    % iterate through stimuli
    for ii = 1 : n_stimuli
        stim_flags = rt_stimuli == stim_set(ii);
        
        % iterate through correctness
        for cc = [0,1]
            correct_flags = (rt_ch == (rt_stimuli > nanmedian(rt_stimuli))) == cc;
            trial_flags = ...
                rt_valid_flags & ...
                contrast_flags & ...
                stim_flags & ...
                correct_flags;
            rts(kk,cc+1).x(ii,1) = normstim_set(ii);
            rts(kk,cc+1).y(ii,1) = nanmedian(rt(trial_flags));
            rts(kk,cc+1).n(ii,1) = sum(trial_flags);
            rts(kk,cc+1).err(ii,:) = ...
                quantile(rt(trial_flags),[.25,.75]) - nanmedian(rt(trial_flags));
        end
    end
end

%% plot reaction times

% figure initialization
fig = figure(figopt,...
    'name',sprintf('reaction_times_%s',contrast_str));

% axes initialization
axes(...
    axesopt.default,...
    axesopt.stimulus);
title('Reaction times');
xlabel(stim_lbl);
ylabel('Reaction time (ms)');

% graphical object preallocation
p = gobjects(n_contrasts,1);

% correctness markers
markers = {'s','o'};

% iterate through contrasts
for kk = 1 : n_contrasts

    % iterate through correctness
    for cc = 1 % [0,1]
        
        % plot error bars
        errorbar(rts(kk,cc+1).x,rts(kk,cc+1).y,...
            rts(kk,cc+1).err(:,1),rts(kk,cc+1).err(:,2),...
            'color',contrast_clrs(kk,:),...
            'marker','none',...
            'linestyle','-',...
            'linewidth',1,...
            'capsize',2);
        
        % plot averages
        p(kk) = plot(rts(kk,cc+1).x,rts(kk,cc+1).y,...
            'color',contrast_clrs(kk,:),...
            'marker',markers{cc+1},...
            'markerfacecolor',contrast_clrs(kk,:),...
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
plot([1,1]*median(normstim_set),ylim,':k');

% legend
leg_str = cellfun(@(x,y)sprintf('%s = %i %s',x,y,contrast_units),...
    repmat({contrast_lbl},n_contrasts,1),num2cell(contrast_set),...
    'uniformoutput',false);
legend(p(isgraphics(p)),leg_str(isgraphics(p)),...
    'position',[0.085,0.66,.27,.2],...
    'box','on');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end