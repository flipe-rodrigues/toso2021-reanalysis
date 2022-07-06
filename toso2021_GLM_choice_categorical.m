%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% design matrix

% number of trials to look back
k = 4;

%
x_i1 = [nan; i1(1:end-1)];
x_i2 = [nan; i2(1:end-1)];
x_t1 = [nan; t1(1:end-1)];
x_t2 = [nan; t2(1:end-1)];
x_ch = [nan; choices(1:end-1)];
x_rwd = [nan; correct(1:end-1)];

%
padded_i1 = [zeros(k - 1,1); x_i1];
padded_i2 = [zeros(k - 1,1); x_i2];
padded_t1 = [zeros(k - 1,1); x_t1];
padded_t2 = [zeros(k - 1,1); x_t2];
padded_c = [zeros(k - 1,1); x_ch];
padded_r = [zeros(k - 1,1); x_rwd];

%
X_i1 = hankel(padded_i1(1 : end - k + 1), x_i1(end - k + 1 : end));
X_i2 = hankel(padded_i2(1 : end - k + 1), x_i2(end - k + 1 : end));
X_t1 = hankel(padded_t1(1 : end - k + 1), x_t1(end - k + 1 : end));
X_t2 = hankel(padded_t2(1 : end - k + 1), x_t2(end - k + 1 : end));
X_c = hankel(padded_c(1 : end - k + 1), x_ch(end - k + 1 : end));
X_r = hankel(padded_r(1 : end - k + 1), x_rwd(end - k + 1 : end));

%
X_i1 = i1 == i_set';
X_i2 = i2 == i_set';
X_t1 = t1 == t_set';
X_t2 = t2 == t_set';

% interactions
t1t2 = t1 .* t2;
t1t2_set = unique(t1t2(valid_flags));
i1i2 = i1 .* i2;
i1i2_set = unique(i1i2(valid_flags));
X_t1t2 = t1t2 == t1t2_set';
X_i1i2 = i1i2 == i1i2_set';
X_t2i2 = t2 .* i2 == unique(t2(valid_flags) .* i2(valid_flags))';
X_t2i1 = t2 .* i1 == unique(t2(valid_flags) .* i1(valid_flags))';
X_t1i2 = t1 .* i2 == unique(t1(valid_flags) .* i2(valid_flags))';
X_t1i1 = t1 .* i1 == unique(t1(valid_flags) .* i1(valid_flags))';

% trial history
X_i1_p = prev_i1 == i_set';
X_i2_p = prev_i2 == i_set';
X_t1_p = prev_t1 == t_set';
X_t2_p = prev_t2 == t_set';

X_c_p = [nan; choices(1:end-1)];
X_c_p2 = [nan; X_c_p(1:end-1)];

X_t1t2_p = [nan;t1t2(1:end-1)] == t1t2_set';
X_i1i2_p = [nan;i1i2(1:end-1)] == i1i2_set';

X_i1_p2 = [nan;prev_i1(1:end-1)] == i_set';
X_i2_p2 = [nan;prev_i2(1:end-1)] == i_set';
X_t1_p2 = [nan;prev_t1(1:end-1)] == t_set';
X_t2_p2 = [nan;prev_t2(1:end-1)] == t_set';

trial_kernel = expkernel('mus',5,'binwidth',1);
i_hist = conv(-prev_i1+prev_i2,trial_kernel.pdf,'same');
t_hist = conv(-prev_t1+prev_t2,trial_kernel.pdf,'same');

%
design_table = table(...X_i1,X_i2,X_t1,X_t2,X_c,X_r,...
    X_i1,X_i1_p,...
    X_i2,X_i2_p,...
    X_t1,X_t1_p,X_t1_p2,...
    X_t2,X_t2_p,X_t2_p2,...
    ...X_i1i2,...
    ...X_t1t2,...,X_t1i1_t,X_t2i2_t...
    ...X_i1_p,X_i2_p,X_t1_p,X_t2_p,...
    ...i_hist,t_hist,...
    ...X_t1t2_p,X_i1i2_p,...
    ...X_i1_p2,X_i2_p2,X_t1_p2,X_t2_p2,...
    X_c,..._p,X_c_p2,...
    X_r...
    );
design = double(design_table.Variables);
coeff_names = ['intercept',design_table.Properties.VariableNames];
n_coeffs = size(design,2) + 1;

%% feature normalization

% z-scoring
mus = nanmean(design,1);
sigs = nanstd(design,0,1);
zdesign = (design - mus) ./ sigs;

%% response variable

% construct response variable
response = choices;

%% trial selection
trial_flags = ...
    ...i1 == i_set(i1_mode_idx) & ...
    valid_flags;
design(~trial_flags,:) = nan;
zdesign(~trial_flags,:) = nan;
response(~trial_flags) = nan;
n_trials = size(design,1);

%% GLM

% distribution selection
distro = 'binomial';

% fit GLM
[B,mdlinfo] = lassoglm(...
    zdesign,response,distro,...
    'standardize',true,...
    'lambda',1e-1,...
    'alpha',1e-2,...
    'CV',10);
[~,nullinfo] = lassoglm(...
    zdesign*0,response,distro,...
    'standardize',true,...
    'lambda',1e-1,...
    'alpha',1e-2,...
    'CV',10);

% extract coefficients~
coeffs = [...
    mdlinfo.Intercept(mdlinfo.IndexMinDeviance);...
    B(:,mdlinfo.IndexMinDeviance)...
    ];

% no regularization
% mdl = fitglm(zdesign,response,...
%     'distribution',distro);
% coeffs = mdl.Coefficients.Estimate;

% compute pseudo r-squared
pseudo_r2 = (nullinfo.Deviance - mdlinfo.Deviance) / nullinfo.Deviance;
pseudo_stimuli = [ones(n_trials,1),zdesign] * coeffs;

%% figure initialization
xxticklabels = [...
    '\beta_0',...
    arrayfun(@(x)sprintf('%s=%i',upper(coeff_names{2}(3:end)),x),i_set','uniformoutput',false),...
    ...arrayfun(@(x)sprintf('%s=%i',upper(coeff_names{3}(3:end)),x),i_set','uniformoutput',false),...
    arrayfun(@(x)sprintf('%s=%i',upper(coeff_names{4}(3:end)),x),i_set','uniformoutput',false),...
    arrayfun(@(x)sprintf('%s=%i',upper(coeff_names{6}(3:end)),x),t_set','uniformoutput',false),...
    arrayfun(@(x)sprintf('%s=%i',upper(coeff_names{9}(3:end)),x),t_set','uniformoutput',false),...
    arrayfun(@(x)sprintf('%s=%i',upper(coeff_names{12}(3:end)),x),i1i2_set','uniformoutput',false),...
    arrayfun(@(x)sprintf('%s=%i',upper(coeff_names{13}(3:end)),x),t1t2_set','uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{end-1}(3:end)),x),k:-1:0+1,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{end}(3:end)),x),k:-1:0+1,'uniformoutput',false),...
    ...arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{2}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    ...arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{3}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    ...arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{4}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    ...arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{5}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    ...arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{6}(3:end)),x),k:-1:0+1,'uniformoutput',false),...
    ...arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{7}(3:end)),x),k:-1:0+1,'uniformoutput',false),...
    ];
figure('name',mfilename,...
    'windowstate','maximized',...
    'numbertitle','off',...
    'inverthardcopy','off',...
    'color','w');
set(gca,axesopt.default,...
    'xlim',[1,n_coeffs]+[-1,1],...
    'xtick',1:n_coeffs,...
    'xticklabel',xxticklabels,...
    'xticklabelrotation',45,...
    'ylim',[0,1]+[-1,1]*.05,...
    'ytick',0,...[0,1],...
    'ticklabelinterpreter','tex',...
    'plotboxaspectratio',[4,1,1],...
    'clipping','off',...
    'ycolor','k');
title(sprintf('C_{t}~%s(\\phi(\\betaX))',distro));
ylabel('$\beta_i$',...
    'interpreter','latex');

%% plot fit coefficients

% plot coefficient relationships
coeff_idcs_offset = 0;
coeff_x_offset = 0;
for jj = 1 : numel(coeff_names)
    coeff_name = coeff_names{jj};
    if jj > 1
        prev_coeff_name = coeff_names{jj-1};
    end
    if strcmpi(coeff_name,'intercept')
        coeff_size = 1;
    else
        coeff_size = size(design_table.(coeff_name),2);
    end
    coeff_idcs = (1 : coeff_size) + coeff_idcs_offset;
    coeff_idcs_offset = coeff_idcs_offset + coeff_size;
    if jj > 1 && contains(coeff_name,prev_coeff_name)
        coeff_x_offset = coeff_x_offset + coeff_size;
%         markersize = markersize * .75;
        color = color + [1,1,1] * .45;
    else
        markersize = 8.5;
        color = [0,0,0];
    end
    coeff_x = coeff_idcs - coeff_x_offset;
    p = plot(coeff_x,coeffs(coeff_idcs),...
        'color',color,...
        'marker','o',...
        'markersize',markersize,...
        'markeredgecolor',color,...
        'markerfacecolor','w',...
        'linewidth',1.5);
    if jj > 1 && contains(coeff_name,prev_coeff_name)
        uistack(p,'bottom');
    end
end

% update axes
axis tight;
xlim(xlim+[-1,1]);
ylim(ylim+[-1,1]*.05*range(ylim));

% zero line
p = plot(xlim,[0,0],'--k',...
    'hittest','off');
uistack(p,'bottom');

% iterate through coefficients
for bb = 1 : max(xlim) - 1
    p = plot([1,1]*bb,ylim,':k',...
        'hittest','off');
    uistack(p,'bottom');
end

% r-squared annotation
text(.95,.05,sprintf('pseudo-R^{2} = %.2f',pseudo_r2),...
    'fontsize',12,...
    'color','k',...
    'horizontalalignment','right',...
    'verticalalignment','bottom',...
    'units','normalized');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end