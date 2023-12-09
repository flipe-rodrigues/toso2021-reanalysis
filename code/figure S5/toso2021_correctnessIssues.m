%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% general settings
m = 1e3;
t = linspace(0,3,m) * max(t_set);
pdf_cutoff = 5 / m;
cdf_cutoff = .01;

%% percept distributions

% preallocation
percept = struct();

% percept definition
percept.mus = t_set';
percept.web = .2;
percept.sig = .05 * max(t_set);
percept.pdfs = normpdf(t,percept.mus',percept.sig);
for tt = 1 : n_t
    percept.pdfs(tt,:) = zeros(1,m);
    [~,dirac_idx] = min(abs(t' - percept.mus(tt)));
    percept.pdfs(tt,dirac_idx) = 1;
    %     percept.pdfs(tt,:) = normpdf(x,percept.mus(tt),percept.sig);
end
percept.pdfs = percept.pdfs ./ nansum(percept.pdfs,2);

% define gaussian kernel to introduce scalar timing
kernel.win = t;
kernel.mus = kernel.win';
kernel.sigs = kernel.mus * percept.web;
kernel.pdfs = normpdf(kernel.win,kernel.mus,kernel.sigs);
I = eye(m);
for ii = 1 : m
    progressreport(ii,m,'generating percept distros');
    if all(isnan(kernel.pdfs(ii,:)))
        kernel.pdfs(ii,:) = I(ii,:);
    end
end
kernel.pdfs = kernel.pdfs ./ nansum(kernel.pdfs,2);

% smear percept distros with gaussian kernel
percept.pdfs = percept.pdfs * kernel.pdfs';
percept.pdfs = percept.pdfs ./ mean(nansum(percept.pdfs,2));
percept.cdfs = cumsum(percept.pdfs,2);

%% color settings
slow_clr = [.0,.4,.95];
avg_clr = [1,1,1] * .0;
fast_clr = [.95,.25,.35];
bg_clr = [1,1,1] * 245 / 255;
clrmap = colorlerp([bg_clr;slow_clr;avg_clr;fast_clr;bg_clr],m);
% clrmap = colorlerp([slow_clr;avg_clr;fast_clr],m);

%% temporal scaling settings

% preallocation
speed = struct();

% speed definition
speed.mus = t';
speed.web = percept.web;
speed.sig = percept.sig;
speed.pdfs = eye(m) * kernel.pdfs';
speed.pdfs = speed.pdfs ./ nansum(speed.pdfs,2);
speed.cdfs = cumsum(speed.pdfs,2);

%% scaling diagram (linear)

% figure initialization
fig = figure(figopt,...
    'color',bg_clr,...
    'name','scaling_marginal_linear');

% axes initialization
xxlim = round([t_set(1),t_set(end)] + ...
    [-1,0] * t_set(1)/range(t_set) * range(t_set));
yylim = xxlim .* [1,2];
xxtick = unique([0,xxlim,t_set']);
yytick = unique([0,yylim,t_set']);
xxticklabel = num2cell(xxtick);
yyticklabel = num2cell(yytick);
axes(axesopt.default,...
    'xlim',xxlim,...
    'ylim',yylim,...
    'xtick',xxtick,...
    'ytick',yytick,...
    'xticklabel',xxticklabel,...
    'yticklabel',yyticklabel,...
    'colormap',clrmap,...
    'clipping','off');
xlabel('Time since stimulus onset (ms)');
ylabel('Internal time since stimulus onset (ms)');

% underlying temporal scaling
x_flags = ...
    t >= xxlim(1) & ...
    t <= xxlim(2);
y_flags = ...
    t >= yylim(1) & ...
    t <= yylim(2);
imagesc(t(x_flags),t(y_flags),speed.cdfs(x_flags,y_flags)');

% categorical boundary
plot([1,1].*t_set(t2_mode_idx),ylim,...
    'color','k',...
    'linestyle',':');
plot(xlim,[1,1].*t_set(t2_mode_idx),...
    'color','k',...
    'linestyle',':');

% iterate through stimuli
for ii = 1 : n_t
    pdf_flags = ...
        (percept.pdfs(ii,:) / max(percept.pdfs(ii,:))) >= pdf_cutoff;
    
    % plot percept distribution
    xpatch = [t(pdf_flags),fliplr(t(pdf_flags))];
    ypatch = [zeros(1,sum(pdf_flags)),fliplr(percept.pdfs(ii,pdf_flags))];
    ypatch = normalize01(ypatch,2) * .1 * max(percept.mus);
    patch(t_set(ii)-ypatch,xpatch,'w',...
        'edgecolor','k',...
        'facealpha',1,...
        'linewidth',1.5,...
        'linestyle','-');
end

% colorbar
clrbar = colorbar();
clrlabel = struct();
clrlabel.string = {'Speed of striatal dynamics (a.u.)'};
clrlabel.fontsize = axesopt.default.fontsize;
clrlabel.rotation = 270;
clrlabel.position = [2.75,.5,0];
set(clrbar,...
    axesopt.colorbar,...
    'ticks',[],...
    'box','on');
set(clrbar.Label,...
    clrlabel);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% scaling diagram (log)

% transfer function
tfun = @(x) log(x);
invfun = @(x) exp(x);

% figure & axes initialization
fig = figure(figopt,...
    'color',bg_clr,...
    'name','scaling_marginal_log');

% axes initialization
xxlim = tfun([t_set(1),t_set(end)]) + [-1,0] * .25 * range(tfun(t_set));
yylim = tfun([t_set(1),t_set(end)]) + [-1,1] * .25 * range(tfun(t_set));
xxtick = unique([xxlim,tfun(t_set)']);
yytick = unique([yylim,tfun(t_set)']);
xxticklabel = num2cell(round(invfun(xxtick)));
yyticklabel = num2cell(round(invfun(yytick)));
axes(axesopt.default,...
    'xlim',xxlim,...
    'ylim',yylim,...
    'xtick',xxtick,...
    'ytick',yytick,...
    'xticklabel',xxticklabel,...
    'yticklabel',yyticklabel,...
    'colormap',clrmap,...
    'clipping','off');
xlabel('Time since stimulus onset (ms)');
ylabel('Internal time since stimulus onset (ms)');

% underlying temporal scaling
x_flags = ...
    t >= invfun(xxlim(1)) & ...
    t <= invfun(xxlim(2));
y_flags = ...
    t >= invfun(yylim(1)) & ...
    t <= invfun(yylim(2));
[T1,T2] = meshgrid(tfun(t(x_flags)),tfun(t(y_flags)));
P = zeros(size(T1));
C = speed.cdfs(x_flags,y_flags)';
surf(T1,T2,P,C,...
    'edgecolor','none');

% categorical boundary
plot(tfun([1,1].*t_set(t2_mode_idx)),ylim,...
    'color','k',...
    'linestyle',':');
plot(xlim,tfun([1,1].*t_set(t2_mode_idx)),...
    'color','k',...
    'linestyle',':');

% iterate through stimuli
for ii = 1 : n_t
    pdf_flags = ...
        (percept.pdfs(ii,:) / max(percept.pdfs(ii,:))) >= pdf_cutoff;
    
    % plot percept distribution
    cdf = cumsum(percept.pdfs(ii,:));
    xpatch = [t(pdf_flags),fliplr(t(pdf_flags))];
    ypatch = [zeros(1,sum(pdf_flags)),fliplr(percept.pdfs(ii,pdf_flags))];
    ypatch = normalize01(ypatch,2) * .05 * tfun(max(percept.mus));
    patch(tfun(t_set(ii))-ypatch,tfun(xpatch),'w',...
        'edgecolor','k',...
        'facealpha',1,...
        'linewidth',1.5,...
        'linestyle','-');
end

% colorbar
clrbar = colorbar();
clrlabel = struct();
clrlabel.string = {'Speed of striatal dynamics (a.u.)'};
clrlabel.fontsize = axesopt.default.fontsize;
clrlabel.rotation = 270;
clrlabel.position = [2.75,.5,0];
set(clrbar,...
    axesopt.colorbar,...
    'ticks',[],...
    'box','on');
set(clrbar.Label,...
    clrlabel);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% joint distribution

% transfer function
tfun = @(x) (x);
invfun = @(x) (x);

% pair specification
t_pairs = [t1,t2];
t_pairset = unique(t_pairs(valid_flags,:),'rows');
n_t_pairs = size(t_pairset,1);

% pair selection
t_pairs2plot = [3,4];

% iterate through S1-S2 pairs
for ii = 1 : n_t_pairs
    t1_idx = find(ismember(t_set,t_pairset(ii,1)));
    t2_idx = find(ismember(t_set,t_pairset(ii,2)));
    
    % figure initialization
    fig = figure(figopt,...
        'name',sprintf('scaling_joint_t1%i_t2%i',t_set(t1_idx),t_set(t2_idx)),...
        'color',bg_clr);
    
    % axes initialization
    xxlim = tfun([0,t_set(end)]) + [-1,1] * tfun(max(t_set)) * .1;
    %     xxlim = tfun([t_set(1),t_set(end)]) + [-1,1] * .1 * range(tfun(t_set));
    yylim = xxlim;
    xxtick = tfun(unique([0,t_set']));
    yytick = xxtick;
    xxticklabel = num2cell(xxtick);
    yyticklabel = xxticklabel;
    xxticklabel(~ismember(xxtick,[0,t_set(end),t_pairset(ii,1)'])) = {''};
    yyticklabel(~ismember(yytick,[0,t_set(end),t_pairset(ii,2)'])) = {''};
    axes(axesopt.default,...
        'xlim',xxlim,...
        'ylim',yylim,...
        'xtick',xxtick,...
        'ytick',yytick,...
        'xticklabel',xxticklabel,...
        'yticklabel',yyticklabel,...
        'clipping','off');
    xlabel(sprintf('%s (%s)',s1_lbl,s_units));
    ylabel(sprintf('%s (%s)',s2_lbl,s_units));
    
    % plot joint distribution
    joint_pdf = ...
        percept.pdfs(t1_idx,:) .* ...
        percept.pdfs(t2_idx,:)';
    joint_pdf = joint_pdf / nansum(joint_pdf,'all');
    [T1,T2] = meshgrid(tfun(t),tfun(t));
    P = normalize01(joint_pdf,[1,2]);
    contourf(T1,T2,P,[1,1]*pdf_cutoff,...
        'color','k',...
        'facecolor','w',...
        'linewidth',1.5);
    
    % iterate through correctness
    for cc = [1,0]
        correct_flags = ...
            P >= pdf_cutoff & ...
            ((invfun(T1) >= invfun(T2)) == (t_pairset(ii,1) >= t_pairset(ii,2))) == cc;
        temp = joint_pdf;
        temp(~correct_flags) = nan;
        xpatch = [tfun(t),fliplr(tfun(t))];
        lim_flags = ...
            xpatch >= xxtick(1) & ...
            xpatch <= xxtick(end);
        
        % plot S1 percept distribution
        s1_marginal_full = nansum(joint_pdf,1);
        s1_marginal = nansum(temp,1);
        s1_marginal = ...
            (s1_marginal - min(s1_marginal_full)) ./ range(s1_marginal_full);
        s1_marginal = normalize01(s1_marginal);
        ypatch = [zeros(size(s1_marginal)),fliplr(s1_marginal)];
        ypatch = ypatch * .1 * tfun(max(percept.mus));
        patch(xpatch(lim_flags),ypatch(lim_flags)+min(ylim),'w',...
            'edgecolor',reward_clrs(cc+1,:),...
            'facealpha',1,...
            'linewidth',1.5,...
            'linestyle','-');
        
        % plot S2 percept distribution
        s2_marginal_full = nansum(joint_pdf,2);
        s2_marginal = nansum(temp,2)';
        s2_marginal = ...
            (s2_marginal - min(s2_marginal_full)) ./ range(s2_marginal_full);
        s2_marginal = normalize01(s2_marginal);
        ypatch = [zeros(size(s2_marginal)),fliplr(s2_marginal)];
        ypatch = ypatch * .1 * tfun(max(percept.mus));
        patch(min(ylim)+ypatch(lim_flags),xpatch(lim_flags),'w',...
            'edgecolor',reward_clrs(cc+1,:),...
            'facealpha',1,...
            'linewidth',1.5,...
            'linestyle','-');
    end
    
    % plot reference lines
    plot(xlim,ylim,':k',...
        'linewidth',1);
    
    % plot NSD lines
    plot(tfun(t_set([1,end-2])),tfun(t_set([1+2,end])),':k',...
        'linewidth',1);
    plot(tfun(t_set([1+2,end])),tfun(t_set([1,end-2])),':k',...
        'linewidth',1);
    
    % plot stimulus pairs
    scatter(...
        tfun(t_pairset(:,1)),...
        tfun(t_pairset(:,2)),...
        50,'w','o','filled',...
        'markeredgecolor',[1,1,1]*.75,...
        'linewidth',1.5);
    scatter(...
        tfun(t_pairset(ii,1)),...
        tfun(t_pairset(ii,2)),...
        50,'w','o','filled',...
        'markeredgecolor','k',...
        'linewidth',1.5);
    
    % annotate NSD = 0
    text(.95,.95,sprintf('%s = 0',nsd_lbl),...
        'units','normalized',...
        'color','k',...
        'rotation',45,...
        'fontsize',axesopt.default.fontsize,...
        'horizontalalignment','right',...
        'verticalalignment','bottom');
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end

%% pseudo decoding

% figure initialization
fig = figure(figopt,...
    'name','pseudo_decoding',...
    'color',bg_clr);

%
clrs = choice_clrs;
clrmap = colorlerp([clrs(1,:);bg_clr;clrs(2,:)],2^8);

% axes initialization
xxlim = round([t_set(1),t_set(end)] + ...
    [-1,0] * t_set(1)/range(t_set) * range(t_set));
yylim = xxlim .* [1,1];
xxtick = unique([0,xxlim,t_set']);
yytick = unique([0,yylim,t_set']);
xxticklabel = num2cell(xxtick);
yyticklabel = num2cell(yytick);
axes(axesopt.default,...
    'xlim',xxlim,...
    'ylim',yylim,...
    'xtick',xxtick,...
    'ytick',yytick,...
    'xticklabel',xxticklabel,...
    'yticklabel',yyticklabel,...
    'colormap',clrmap,...
    'clipping','off');
xlabel(sprintf('Time since %s onset (%s)',s2_lbl,s_units));
ylabel(sprintf('Internal time since %s onset (%s)',s2_lbl,s_units));

%
[T1,T2] = meshgrid(t,t);

% preallocation
p_tR = zeros(m,m,n_choices);

% iterate through S1-S2 pairs
for ii = 1 : n_t_pairs
    t1_idx = find(ismember(t_set,t_pairset(ii,1)));
    t2_idx = find(ismember(t_set,t_pairset(ii,2)));
    t1_flags = t <= t_set(t1_idx);
    t2_flags = t <= t_set(t2_idx);

% % iterate through S2 durations
% for ii = 1 : n_t
%     
    % iterate through correctness
    for cc = [1,0]
        choice_flags = ...
            ((T1 >= T2) == cc);
        correct_flags = ...
            ((T1 >= T2) == (t_pairset(ii,1) >= t_pairset(ii,2))) == cc;
        temp = speed.pdfs .* choice_flags;
        temp(~t2_flags,:) = 0;
        p_tR(:,:,cc+1) = p_tR(:,:,cc+1) + temp / n_t_pairs;
    end
end

% underlying temporal scaling
x_flags = ...
    t >= xxlim(1) & ...
    t <= xxlim(2);
y_flags = ...
    t >= yylim(1) & ...
    t <= yylim(2);
p_diff = diff(p_tR(x_flags,y_flags,:),1,3);
p_1 = p_tR(x_flags,y_flags,1);
p_2 = p_tR(x_flags,y_flags,2);
% p_1 = p_1 ./ nansum(p_1,2);
% p_2 = p_2 ./ nansum(p_2,2);
% p_diff = p_2 - p_1;
imagesc(t(x_flags),t(y_flags),p_diff',[-1,1]*25/m);

% convert from tensor to rgb
P_tR_avg = p_tR(x_flags,y_flags,:);
P_tR_avg = P_tR_avg ./ nansum(P_tR_avg,2);
P = tensor2rgb(permute(P_tR_avg,[2,1,3]),choice_clrs);
% imagesc(x(x_flags),x(y_flags),P);

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end