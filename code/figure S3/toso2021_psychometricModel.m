%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% reversed cumulative gumbel definition
z = @(x) log(-log(x));
sigmoidrgumbel = @(x, m, w, a) ...
    exp(-exp((z(1 - a) - z(a)) / w * (x - m) + z(.5)));

%% psychometric function definition
psi = @(x, m, w, a, g, l, S) ...
    g + (1 - l - g) * S(x,m,w,a);

%% plot psychometric model specification

% figure & axes initialization
fig = figure(figopt,...
    'name','psychometric_model');

% axes settings
axes(...
    axesopt.default,...
    axesopt.stimulus,...
    axesopt.psycurve);

% labels
xlabel(sprintf('%s (%s)',s2_lbl,s_units));
ylabel(sprintf('P(%s > %s)',s2_lbl,s1_lbl));

% stimulus level
x = linspace(-.05,1.05,1e3);

% psychometric function parameters
m = .3;
w = .8;
a = .05;
g = .1;
l = .15;

% compute sigmoid
S = sigmoidrgumbel(x,m,w,a);

% compute psychometric curve
y = psi(x,m,w,a,g,l,sigmoidrgumbel);

% plot psychometric curve
plot(x,y,...
    'color','k',...
    'linewidth',1.5);

% plot unscaled sigmoid
plot(x,S,...
    'color','k',...
    'linestyle','--',...
    'linewidth',1);

% plot threshold
plot([1,1]*m,[0,.5],...
    'color','k');
plot([min(xlim),m],[1,1]*.5,...
    'color','k',...
    'linestyle',':');

% plot width
[~,w1_idx] = min(abs(S - a));
[~,w2_idx] = min(abs(S - (1-a)));
w1 = round(x(w1_idx),2);
w2 = round(x(w2_idx),2);
plot([w1,w2],[0,0],...
    'color','k');
plot([1,1]*w1,[0,1]*a,...
    'color','k');
plot([1,1]*w2,[0,1]*(1-a),...
    'color','k');
plot([min(xlim),w1],[1,1]*a,...
    'color','k',...
    'linestyle',':');
plot([min(xlim),w2],[1,1]*(1-a),...
    'color','k',...
    'linestyle',':');

% plot lapse rate
plot(xlim,[1,1]*(1-l),...
    'color','k',...
    'linestyle',':');

% legend
leg = legend('$\psi(x;m,w(\alpha),\gamma,\lambda)$',...
    '$S(x;m,w(\alpha))$');
set(leg,...
    'box','off',...
    'units','normalized',...
    'position',[.235,.7,leg.Position(3:end)],...
    'fontsize',12,...
    'linewidth',1,...
    'interpreter','latex');

% update axes
xxtick = unique([normstim_set',m]);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,[normstim_set([1,end])',m])) = {''};
xxticklabel(xxtick == m) = {'$m$'};
% xxticklabel(xxtick == 1) = {'$w(\alpha)$'};
yytick = unique([0,.5,1,a,1-a,g,1-l]);
yyticklabel = num2cell(yytick);
yyticklabel(yytick == a) = {'$\alpha$'};
yyticklabel(yytick == 1 - a) = {'$1-\alpha$'};
yyticklabel(yytick == g) = {'$\gamma$'};
yyticklabel(yytick == 1 - l) = {'$1-\lambda$'};
set(gca,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ytick',yytick,...
    'yticklabel',yyticklabel,...
    'ticklabelinterpreter','latex');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end