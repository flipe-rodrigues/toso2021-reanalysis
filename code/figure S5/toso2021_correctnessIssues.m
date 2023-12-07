%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% general settings
m = 1e3;
x = linspace(-1,5,m) * max(t_set);
pdf_cutoff = 1e-5;

%% percept distributions

% define percept distributions
percept = struct();

% percept definition
percept.mus = t_set';
percept.web = .2;
percept.sig = 0;
percept.pdfs = normpdf(x,percept.mus',percept.sig);
for ii = 1 : n_t
    percept.pdfs(ii,:) = zeros(1,m);
    [~,dirac_idx] = min(abs(x' - percept.mus(ii)));
    percept.pdfs(ii,dirac_idx) = 1;
end
percept.pdfs = percept.pdfs ./ nansum(percept.pdfs,2);

% define gaussian kernel to introduce scalar timing
kernel.win = x;
kernel.mus = kernel.win';
kernel.sigs = kernel.mus * percept.web;
kernel.pdfs = normpdf(kernel.win,kernel.mus,kernel.sigs);
I = eye(m);
for jj = 1 : m
    if all(isnan(kernel.pdfs(jj,:)))
        kernel.pdfs(jj,:) = I(jj,:);
    end
end
kernel.pdfs = kernel.pdfs ./ nansum(kernel.pdfs,2);

% smear percept distros with gaussian kernel
percept.pdfs = percept.pdfs * kernel.pdfs';
percept.pdfs = percept.pdfs ./ nansum(percept.pdfs,2);
percept.cdfs = cumsum(percept.pdfs,2);

%% color settings
slow_clr = [.0,.4,.95];
avg_clr = [1,1,1] * .0;
fast_clr = [.95,.25,.35];
clrmap = colorlerp([slow_clr;avg_clr;fast_clr],m);

%% temporal scaling settings
speed.mus = x';
speed.web = percept.web;
speed.sig = percept.sig;
speed.pdfs = normpdf(x,speed.mus',speed.sig);
for ii = 1 : m
    speed.pdfs(ii,:) = zeros(1,m);
    [~,dirac_idx] = min(abs(x' - speed.mus(ii)));
    speed.pdfs(ii,dirac_idx) = 1;
end
speed.pdfs = speed.pdfs ./ nansum(speed.pdfs,2);
speed.pdfs = speed.pdfs * kernel.pdfs';
speed.pdfs = speed.pdfs ./ nansum(speed.pdfs,2);
speed.cdfs = cumsum(speed.pdfs,2);

%%
x_alphas = linspace(0,2,m);
y_alphas = zeros(1,m);
[~,dirac_idx] = min(abs(x_alphas - 1));
y_alphas(dirac_idx) = 1;
y_alphas = y_alphas * kernel.pdfs'; % normalize01(normpdf(linspace(-1,1,m),0,percept.web*2),2);
y_alphas = normalize01(y_alphas);

slopes = nan(m,1);
for ii = 1 : m
    idx = find(rand >= slope_cdf,1);
end
idcs = sum(rand(m,1) >= slope_cdf',2);
slopes = sort(x_alphas(idcs));

figure; plot(x_alphas,y_alphas,x_alphas,slope_cdf);
%% scaling diagram (linear)

% figure initialization
fig = figure(figopt,...
    'name','scaling_diagram_linear');

% axes initialization
xxlim = round([t_set(1),t_set(end)] + ...
    [-1,0] * t_set(1)/range(t_set) * range(t_set));
axes(axesopt.default,...
    'xlim',xxlim,...
    'ylim',xxlim.*[1,1.5],...
    'xtick',[0;t_set],...
    'ytick',[0;t_set],...
    'xticklabel',['0';axesopt.stimulus.xticklabel],...
    'yticklabel',['0';axesopt.stimulus.xticklabel],...
    'colormap',clrmap,...
    'clipping','on');
xlabel('Time since stimulus onset (s)');
ylabel('Internal time since stimulus onset (s)');

% underlying temporal scaling
% speed_bounds = [1,1] + [-1,1] * percept.web * 1.05 * 5;
% slopes = linspace(speed_bounds(1),speed_bounds(2),m);
% slopes = 1;
intercepts = linspace(-percept.sig,percept.sig,m) * 5;
xx_bounds = [0,t_set(end)];
xx = linspace(xx_bounds(1),xx_bounds(2),m);
for ii = 1 : m
    yy = xx * slopes(ii) + intercepts(ii);
    plot(xx,yy,...
        'color',[clrmap(ii,:)],...
        'linewidth',.1);
end

% categorical boundary
plot([1,1].*t_set(t2_mode_idx),ylim,...
    'color','k',...
    'linestyle',':');
plot(xlim,[1,1].*t_set(t2_mode_idx),...
    'color','k',...
    'linestyle',':');

% iterate through stimuli
for ii = 1 : n_t
    xflags = percept.pdfs(ii,:) > pdf_cutoff;
    
    % plot percept distribution
    cdf = cumsum(percept.pdfs(ii,:));
    xpatch = [x(xflags),fliplr(x(xflags))];
    ypatch = [zeros(1,sum(xflags)),fliplr(percept.pdfs(ii,xflags))];
    ypatch = normalize01(ypatch,2) * .1 * max(percept.mus);
    if ismember(ii,n_t/2+[0,1])
        edgecolor = 'k';
        facealpha = 1;
        cpatch = [1,1,1];
    else
        edgecolor = 'k';
        facealpha = 1;
        cpatch = [1,1,1];
    end
    patch(t_set(ii)-ypatch,xpatch,cpatch,...
        'edgecolor',edgecolor,...
        'facealpha',facealpha,...
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
    'name','scaling_diagram_log');

% axes initialization
xxlim = tfun([t_set(1),t_set(end)]) + ...
    [-1,0] * t_set(1)/range(t_set) * tfun(range(t_set));
axes(axesopt.default,...
    'xlim',xxlim,...
    'ylimspec','tight',...
    'ylimspec','tight',...
    'xtick',tfun(t_set),...
    'ytick',tfun(t_set),...
    'xticklabel',axesopt.stimulus.xticklabel,...
    'yticklabel',axesopt.stimulus.xticklabel,...
    'colormap',clrmap,...
    'clipping','on');
xlabel('Time since stimulus onset (s)');
ylabel('Internal time since stimulus onset (s)');

% underlying temporal scaling
speed_bounds = [1,1] + [-1,1] * percept.web * 1.05 * 2;
slopes = linspace(speed_bounds(1),speed_bounds(2),m);
intercepts = linspace(-percept.sig,percept.sig,m) * 2;
xx_bounds = [t_set(1),t_set(end)] + [-1,0] * .1 * range(t_set);
xx = linspace(xx_bounds(1),xx_bounds(2),m);
for ii = 1 : m
    yy = xx * slopes(ii) + intercepts(ii);
    plot(tfun(xx),tfun(yy),...
        'color',[clrmap(ii,:),alphas(ii)],...
        'linewidth',.1);
end

% categorical boundary
% plot([1,1].*t_set(t2_mode_idx),ylim,...
%     'color','k',...
%     'linestyle',':');
% plot(xlim,[1,1].*t_set(t2_mode_idx),...
%     'color','k',...
%     'linestyle',':');

% iterate through stimuli
for ii = 1 : n_t
    xflags = percept.pdfs(ii,:) > pdf_cutoff;
    
    % plot percept distribution
    cdf = cumsum(percept.pdfs(ii,:));
    xpatch = [x(xflags),fliplr(x(xflags))];
    ypatch = [zeros(1,sum(xflags)),fliplr(percept.pdfs(ii,xflags))];
    ypatch = normalize01(ypatch,2) * .1 * tfun(max(percept.mus));
    if ismember(ii,n_t/2+[0,1])
        edgecolor = 'k';
        facealpha = 1;
        cpatch = [1,1,1];
    else
        edgecolor = 'k';
        facealpha = 1;
        cpatch = [1,1,1];
    end
    patch(tfun(t_set(ii))-ypatch,tfun(xpatch),cpatch,...
        'edgecolor',edgecolor,...
        'facealpha',facealpha,...
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
s_pairs = [s1,s2];
s_pairset = unique(s_pairs(valid_flags,:),'rows');
n_s_pairs = size(s_pairset,1);

% preallocation
p_choice = nan(n_s_pairs,1);
n_trials_perpair = nan(n_s_pairs,1);

% iterate through S1-S2 pairs
for ii = 1 : n_s_pairs
    s_flags = all(s_pairs == s_pairset(ii,:),2);
    trial_flags = ...
        valid_flags & ...
        unique_flags & ...
        s_flags;
    if sum(trial_flags) == 0
        continue;
    end
    
    % compute average performance for the current pair
    p_choice(ii) = mean(choice(trial_flags));
    n_trials_perpair(ii) = sum(trial_flags);
end

% nan filtering
s_pairset = s_pairset(~isnan(p_choice),:);
p_choice = p_choice(~isnan(p_choice));

% figure & axes initialization
fig = figure(figopt,...
    'name','generalization_matrix_Si');
axes(axesopt.default,...
    'xlim',tfun([s_set(1),s_set(end)]) + [-1,1] * .1 * range(tfun(s_set)),...
    'ylim',tfun([s_set(1),s_set(end)]) + [-1,1] * .1 * range(tfun(s_set)),...
    'xtick',tfun(s_set),...
    'ytick',tfun(s_set),...
    'xticklabel',num2cell(s_set),...
    'yticklabel',num2cell(s_set),...
    'clim',clims,...
    'clipping','off');
xlabel(sprintf('%s (%s)',s1_lbl,s_units));
ylabel(sprintf('%s (%s)',s2_lbl,s_units));

% colorbar
cbar = colorbar(...
    'color','k',...
    'limits',clims,...
    'box',axesopt.default.box,...
    'linewidth',axesopt.default.linewidth,...
    'tickdirection','out',...
    'ticklength',unique(axesopt.default.ticklength),...
    'fontsize',axesopt.default.fontsize,...
    'ticks',unique([clims,.5]));
cbar.Label.String = sprintf('P(%s > %s)',s2_lbl,s1_lbl);
cbar.Label.Rotation = -90;
cbar.Label.VerticalAlignment = 'bottom';

% plot reference lines
plot(xlim,ylim,':k',...
    'linewidth',1);

% plot NSD lines
plot(tfun(s_set([1,end-2])),tfun(s_set([1+2,end])),':k',...
    'linewidth',1);
plot(tfun(s_set([1+2,end])),tfun(s_set([1,end-2])),':k',...
    'linewidth',1);

% iterate through S1-S2 pairs
for ii = [5,6] % 1 : n_s_pairs
    s1_flags = ismember(t_set,s_pairset(ii,1));
    s2_flags = ismember(t_set,s_pairset(ii,2));
    joint_pdf = percept.pdfs(s1_flags,:) .* percept.pdfs(s2_flags,:)';
    
    [X,Y] = meshgrid(tfun(x),tfun(x));
    Z = percept.pdfs(s1_flags,:) .* percept.pdfs(s2_flags,:)';
    Z = Z / nansum(Z,'all');
    contour(X,Y,Z,[1,1]*pdf_cutoff,...
        'color','k',...
        'linewidth',1.5);
end

% plot performance
scatter(...
    tfun(s_pairset(:,1)),...
    tfun(s_pairset(:,2)),...
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