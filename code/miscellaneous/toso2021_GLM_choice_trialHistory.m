%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% design matrix

% number of trials to look back
k = 5;

%
x_i1 = i1;
x_i2 = i2;
x_t1 = t1;
x_t2 = t2;
x_ch = [nan; choice(1:end-1)];
x_rwd = [nan; correct(1:end-1)];
x_i1i2 = i1 .* i2;
x_i1t1 = i1 .* t1;
x_i2t2 = i2 .* t2;
x_t1t2 = t1 .* t2;

%
padded_i1 = [zeros(k - 1,1); x_i1];
padded_i2 = [zeros(k - 1,1); x_i2];
padded_t1 = [zeros(k - 1,1); x_t1];
padded_t2 = [zeros(k - 1,1); x_t2];
padded_c = [zeros(k - 1,1); x_ch];
padded_r = [zeros(k - 1,1); x_rwd];
padded_i1i2 = [zeros(k - 1,1); x_i1i2];
padded_i1t1 = [zeros(k - 1,1); x_i1t1];
padded_i2t2 = [zeros(k - 1,1); x_i2t2];
padded_t1t2 = [zeros(k - 1,1); x_t1t2];

%
X_i1 = hankel(padded_i1(1 : end - k + 1), x_i1(end - k + 1 : end));
X_i2 = hankel(padded_i2(1 : end - k + 1), x_i2(end - k + 1 : end));
X_t1 = hankel(padded_t1(1 : end - k + 1), x_t1(end - k + 1 : end));
X_t2 = hankel(padded_t2(1 : end - k + 1), x_t2(end - k + 1 : end));
X_c = hankel(padded_c(1 : end - k + 1), x_ch(end - k + 1 : end));
X_r = hankel(padded_r(1 : end - k + 1), x_rwd(end - k + 1 : end));
X_i1i2 = hankel(padded_i1i2(1 : end - k + 1), x_i1i2(end - k + 1 : end));
X_i1t1 = hankel(padded_i1t1(1 : end - k + 1), x_i1t1(end - k + 1 : end));
X_i2t2 = hankel(padded_i2t2(1 : end - k + 1), x_i2t2(end - k + 1 : end));
X_t1t2 = hankel(padded_t1t2(1 : end - k + 1), x_t1t2(end - k + 1 : end));

%
design_table = table(X_i1,X_i2,X_t1,X_t2,X_c,X_r,X_i1i2,X_t1t2);
design = design_table.Variables;
coeff_names = ['intercept',design_table.Properties.VariableNames];
n_coeffs = size(design,2) + 1;

%% feature normalization

% z-scoring
mus = nanmean(design,1);
sigs = nanstd(design,0,1);
zdesign = (design - mus) ./ sigs;

%% figure initialization
xxticklabels = [...
    '\beta_0',...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{2}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{3}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{4}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{5}(3:end)),x),k-1:-1:0,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{6}(3:end)),x),k:-1:0+1,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',upper(coeff_names{7}(3:end)),x),k:-1:0+1,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',insertAfter(upper(coeff_names{8}(3:end)),2,':'),x),k-1:-1:0,'uniformoutput',false),...
    arrayfun(@(x)sprintf('%s_{t-%i}',insertAfter(upper(coeff_names{9}(3:end)),2,':'),x),k-1:-1:0,'uniformoutput',false),...
    ];
figure('name',mfilename,...
    'windowstate','maximized',...
    'numbertitle','off',...
    'inverthardcopy','off',...
    'color','w');
pbar = 4;
set(gca,axesopt.default,...
    'xlim',[1,n_coeffs]+[-1,1],...
    'xtick',1:n_coeffs,...
    'xticklabel',xxticklabels,...
    'xticklabelrotation',90,...
    'ylim',[0,1]+[-1,1]*.05,...
    'ytick',0,...[0,1],...
    'ticklabelinterpreter','tex',...
    'plotboxaspectratio',[pbar,1,1],...
    'ticklength',axesopt.default.ticklength/pbar,...
    'clipping','off',...
    'ycolor','k');
title(sprintf('C_{t}~%s(\\phi(\\betaX))',distro));
ylabel('$\beta_i$',...
    'interpreter','latex');

%% response variable

% construct response variable
response = choice;

% iterate through subjects
for ss = 1 : n_subjects + 1
    if ss <= n_subjects
        subject_flags = subjects == subject_set(ss);
    else
        subject_flags = ismember(subjects,subject_set);
    end
    
    % trial selection
    trial_flags = ...
        unique_flags & ...
        data.Trial >= k & ...
        valid_flags & ...
        subject_flags;
%     design(~trial_flags,:) = nan;
%     zdesign(~trial_flags,:) = nan;
%     response(~trial_flags) = nan;
%     n_trials = size(design,1);
    n_trials = sum(trial_flags);
    
    %% GLM
    
    % distribution selection
    distro = 'binomial';
    
    % fit GLM
    lambda = 1e-1;
    alpha = 1e-2;
    [B,mdlinfo] = lassoglm(...
        zdesign(trial_flags,:),response(trial_flags),distro,...
        'standardize',true,...
        'lambda',lambda,...
        'alpha',alpha,...
        'CV',10);
    [~,nullinfo] = lassoglm(...
        zdesign(trial_flags,:)*0,response(trial_flags),distro,...
        'standardize',true,...
        'lambda',lambda,...
        'alpha',alpha,...
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
    pseudo_stimuli = [ones(n_trials,1),zdesign(trial_flags,:)] * coeffs;
    
    %% plot fit coefficients
    
    % plot coefficient relationships
    coeff_offset = 0;
    for jj = 1 : numel(coeff_names)
        coeff_name = coeff_names{jj};
        if strcmpi(coeff_name,'intercept')
            coeff_size = 1;
        else
            coeff_size = size(design_table.(coeff_name),2);
        end
        coeff_idcs = (1 : coeff_size) + coeff_offset;
        coeff_offset = coeff_offset + coeff_size;
        plot(coeff_idcs,coeffs(coeff_idcs),...
            'color',subject_clr*(ss<=n_subjects),...
            'marker','none',...
            'hittest','off',...
            'linewidth',2);
    end
    
    % iterate through coefficients
    for bb = 1 : n_coeffs
        plot(bb,nanmean(coeffs(bb)),...
            'marker','o',...
            'markersize',8.5,...
            'markeredgecolor',subject_clr*(ss<=n_subjects),...
            'markerfacecolor','w',...
            'linewidth',1.5);
    end
    
    % update axes
    axis tight;
    xlim([1,n_coeffs]+[-1,1]);
    ylim(ylim+[-1,1]*.05*range(ylim));
    
    % zero line
    p = plot(xlim,[0,0],'--k',...
        'hittest','off');
    uistack(p,'bottom');
    
    % iterate through coefficients
    for bb = 1 : n_coeffs
        p = plot([1,1]*bb,ylim,':k',...
            'hittest','off');
        uistack(p,'bottom');
    end
    
    % r-squared annotation
    text(.95,.8+.05*ss,sprintf('pseudo-R^{2} = %.2f',pseudo_r2),...
        'fontsize',12,...
        'color',subject_clr*(ss<=n_subjects),...
        'horizontalalignment','right',...
        'verticalalignment','top',...
        'units','normalized');
end