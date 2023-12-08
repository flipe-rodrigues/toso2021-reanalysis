%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% general settings
m = 1e3;
x = linspace(0,3,m) * max(t_set);
pdf_cutoff = 1 / m;
cdf_cutoff = .01;

%% percept distributions

% preallocation
percept = struct();

% percept definition
percept.mus = t_set';
percept.web = .2;
percept.sig = 0;
percept.pdfs = normpdf(x,percept.mus',percept.sig);
for tt = 1 : n_t
    percept.pdfs(tt,:) = zeros(1,m);
    [~,dirac_idx] = min(abs(x' - percept.mus(tt)));
    percept.pdfs(tt,dirac_idx) = 1;
end
percept.pdfs = percept.pdfs ./ nansum(percept.pdfs,2);

% define gaussian kernel to introduce scalar timing
kernel.win = x;
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
percept.pdfs = percept.pdfs ./ nansum(percept.pdfs,2);
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
speed.mus = x';
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
xlabel('Time since stimulus onset (s)');
ylabel('Internal time since stimulus onset (s)');

% underlying temporal scaling
x_flags = ...
    x >= xxlim(1) & ...
    x <= xxlim(2);
y_flags = ...
    x >= yylim(1) & ...
    x <= yylim(2);
imagesc(x(x_flags),x(y_flags),speed.cdfs(x_flags,y_flags)');

% categorical boundary
plot([1,1].*t_set(t2_mode_idx),ylim,...
    'color','k',...
    'linestyle',':');
plot(xlim,[1,1].*t_set(t2_mode_idx),...
    'color','k',...
    'linestyle',':');

% iterate through stimuli
for ii = 1 : n_t
    cdf_flags = ...
        percept.cdfs(ii,:) >= cdf_cutoff & ...
        percept.cdfs(ii,:) <= (1 - cdf_cutoff);
    
    % plot percept distribution
    xpatch = [x(cdf_flags),fliplr(x(cdf_flags))];
    ypatch = [zeros(1,sum(cdf_flags)),fliplr(percept.pdfs(ii,cdf_flags))];
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
xlabel('Time since stimulus onset (s)');
ylabel('Internal time since stimulus onset (s)');

% underlying temporal scaling
x_flags = ...
    x >= invfun(xxlim(1)) & ...
    x <= invfun(xxlim(2));
y_flags = ...
    x >= invfun(yylim(1)) & ...
    x <= invfun(yylim(2));
[X,Y] = meshgrid(tfun(x(x_flags)),tfun(x(y_flags)));
Z = zeros(size(X));
C = speed.cdfs(x_flags,y_flags)';
surf(X,Y,Z,C,...
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
    cdf_flags = ...
        percept.cdfs(ii,:) >= cdf_cutoff & ...
        percept.cdfs(ii,:) <= (1 - cdf_cutoff);
    
    % plot percept distribution
    cdf = cumsum(percept.pdfs(ii,:));
    xpatch = [x(cdf_flags),fliplr(x(cdf_flags))];
    ypatch = [zeros(1,sum(cdf_flags)),fliplr(percept.pdfs(ii,cdf_flags))];
    ypatch = normalize01(ypatch,2) * .05 * tfun(max(percept.mus));
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
t_pairs = [t1,t2];
t_pairset = unique(t_pairs(valid_flags,:),'rows');

% figure & axes initialization
fig = figure(figopt,...
    'color',bg_clr,...
    'position',[744 630 460 470],...
    'name','scaling_joint_log');
axes(axesopt.default,...
    'xlim',tfun([t_set(1),t_set(end)]) + [-1,1] * .1 * range(tfun(t_set)),...
    'ylim',tfun([t_set(1),t_set(end)]) + [-1,1] * .1 * range(tfun(t_set)),...
    'xtick',tfun(t_set),...
    'ytick',tfun(t_set),...
    'xticklabel',num2cell(t_set),...
    'yticklabel',num2cell(t_set),...
    'clim',clims,...
    'clipping','off');
xlabel(sprintf('%s (%s)',s1_lbl,s_units));
ylabel(sprintf('%s (%s)',s2_lbl,s_units));

% plot reference lines
plot(xlim,ylim,':k',...
    'linewidth',1);

% plot NSD lines
plot(tfun(t_set([1,end-2])),tfun(t_set([1+2,end])),':k',...
    'linewidth',1);
plot(tfun(t_set([1+2,end])),tfun(t_set([1,end-2])),':k',...
    'linewidth',1);

% iterate through S1-S2 pairs
for ii = [5,6] % 1 : n_s_pairs
    s1_flags = ismember(t_set,t_pairset(ii,1));
    s2_flags = ismember(t_set,t_pairset(ii,2));
    joint_pdf = percept.pdfs(s1_flags,:) .* percept.pdfs(s2_flags,:)';
    joint_pdf = joint_pdf / nansum(joint_pdf,'all');
    [X,Y] = meshgrid(tfun(x),tfun(x));
    Z = joint_pdf;
    contour(X,Y,Z,[1,1]*pdf_cutoff^2,...
        'color','k',...
        'linewidth',1.5);
end

% plot performance
scatter(...
    tfun(t_pairset(:,1)),...
    tfun(t_pairset(:,2)),...
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

%% example trials (linear)

% figure initialization
fig = figure(figopt,...
    'color',bg_clr,...
    'name','scaling_examples_linear');

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
xlabel('Time since stimulus onset (s)');
ylabel('Internal time since stimulus onset (s)');

% underlying temporal scaling
x_flags = ...
    x >= xxlim(1) & ...
    x <= xxlim(2);
y_flags = ...
    x >= yylim(1) & ...
    x <= yylim(2);
imagesc(x(x_flags),x(y_flags),speed.cdfs(x_flags,y_flags)');

% categorical boundary
plot([1,1].*t_set(t2_mode_idx),ylim,...
    'color','k',...
    'linestyle',':');
plot(xlim,[1,1].*t_set(t2_mode_idx),...
    'color','k',...
    'linestyle',':');

% iterate through stimuli
for ii = 1 : n_t
    cdf_flags = percept.pdfs(ii,:) > pdf_cutoff;
    
    % plot percept distribution
    xpatch = [x(cdf_flags),fliplr(x(cdf_flags))];
    ypatch = [zeros(1,sum(cdf_flags)),fliplr(percept.pdfs(ii,cdf_flags))];
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