%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% simulation settings
N = 20;
T = n_tbins;
S = 2;
K = 100;

%% ROI & time settings
si_roi = [0,t_set(end)];
ti = si_roi(1);
tf = si_roi(2);
t = linspace(ti,tf*2,T);

%% model parameters
mus = [linspace(ti,tf,N/2),linspace(ti,tf,N/2)];
sigma = .05 * range(si_roi);
x = normpdf(t',mus,sigma);

%% temporal scale settings
tsi = 0;
tsf = 2;
ts_support = linspace(tsi,tsf,1e2);
ts_mu = 1;
ts_sigma = .05;
ts_pd = truncate(makedist('normal',...
    'mu',ts_mu,...
    'sigma',ts_sigma),0,inf);
ts_pdf = ts_pd.pdf(ts_support);
figure; plot(ts_support,ts_pdf);

%% general settings
m = 10e3;
x = linspace(-1,5,m);

%% percept distributions

% define percept distributions
percept = struct();

% percept definition
percept.mus = t_set' / max(t_set);
percept.web = .1;
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
clrmap_res = 1e3;
clrmap = colorlerp([slow_clr;avg_clr;fast_clr],clrmap_res);

% figure initialization
fig = figure(figopt,...
    'name','scaling_diagram');
axes(axesopt.default,...
    'xlim',[0,1],...
    'ylim',[0,1]+[-0,1]*.3,...
    'xtick',[0;axesopt.stimulus.xtick],...
    'ytick',[0;axesopt.stimulus.xtick],...
    'xticklabel',['0';axesopt.stimulus.xticklabel],...
    'yticklabel',['0';axesopt.stimulus.xticklabel],...
    'colormap',clrmap,...
    'clipping','off');
xlabel('Time since stimulus onset (s)');
ylabel('Internal time since stimulus onset (s)');

%
speed_bounds = [1,1] + [-1,1] * percept.web * 1.05 * 2;
slopes = linspace(speed_bounds(1),speed_bounds(2),clrmap_res);
intercepts = linspace(-percept.sig,percept.sig,clrmap_res) * 2;
xx = linspace(0,1,clrmap_res);
for ii = 1 : clrmap_res
    yy = xx * slopes(ii) + intercepts(ii);
    clrs = colorlerp([clrmap(ii,:); [1,1,1]],3);
    clr = clrmap(ii,:);
    plot(xx,yy,...
        'color',clrs(end-1,:),...
        'linewidth',.1);
end

% categorical boundary
plot([1,1].*boundary,ylim,...
    'color','k',...
    'linestyle',':');
plot(xlim,[1,1].*boundary,...
    'color','k',...
    'linestyle',':');

% iterate through stimuli
for ii = 1 : n_t
    xflags = percept.pdfs(ii,:) > 1e-5;

    % plot percept distribution
    cdf = cumsum(percept.pdfs(ii,:));
    xpatch = [x(xflags),fliplr(x(xflags))];
    ypatch = [zeros(1,sum(xflags)),fliplr(percept.pdfs(ii,xflags))]*3;
    ypatch = normalize01(ypatch,2) * .05;
    if ismember(ii,n_t/2+[0,1])
        edgecolor = 'k';
        facealpha = 1;
        cpatch = [1,1,1];
    else
        edgecolor = 'k';
        facealpha = 1;
        cpatch = [1,1,1];
    end
    patch(t_set(ii)/max(t_set)-ypatch,xpatch,cpatch,...
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