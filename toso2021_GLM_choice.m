%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% choice GLM (complete model)

% design matrix
X = [s1,s2,d1,d2];
n_regressors = size(X,2);
n_coefficients = n_regressors + 1;

% feature normalization
Z = (X - nanmean(X)) ./ nanstd(X);

% preallocation
betas = nan(n_subjects,n_coefficients);

% iterate through subjects
for ss = 1 : n_subjects
    subject_flags = subjects == subject_set(ss);
    trial_flags = ...
        valid_flags & ...
        subject_flags;
    
    % fit GLM to each subject
    mdl = fitglm(Z(trial_flags,:),choice(trial_flags,:),'linear',...
        'predictorvars',{s1_lbl,s2_lbl,d1_lbl,d2_lbl},...
        'distribution','binomial',...
        'intercept',true);
    betas(ss,:) = mdl.Coefficients.Estimate;
end

% fit GLM to subject pool
mdl = fitglm(Z(valid_flags,:),choice(valid_flags,:),'linear',...
    'predictorvars',{s1_lbl,s2_lbl,d1_lbl,d2_lbl},...,'T1:I1','T2:I2'},...
    'distribution','binomial',...
    'intercept',true);
bigbetas = mdl.Coefficients.Estimate;
beta_s1 = betas(2);
beta_s2 = bigbetas(3);
beta_labels = mdl.CoefficientNames;
beta_labels{1} = 'Intercept';

%% plot GLM coefficients
fig = figure(figopt,...
    'name','choice_GLM');
axes(axesopt.default,...
    'xlim',[0,n_coefficients+1],...
    'xtick',1:n_coefficients,...
    'xticklabelrotation',0,...
    'xticklabel',beta_labels,...
    'plotboxaspectratio',[1,1,1]);
title(sprintf('%s>%s~Binomial(\\phi(\\betaX))',s2_lbl,s1_lbl));
xlabel('Regressor');
ylabel('Weight');

% iterate through subjects
for ss = 1 : n_subjects
   offset = (ss - (n_subjects + 1) / 2) * .025 * range(xlim);
   
    % plot subject coefficients
    stem((1:n_coefficients)+offset,betas(ss,:),...
        'color',subject_clr,...
        'marker','none',...
        'linewidth',1.5);
end

% iterate through coefficients
for bb = 1 : n_coefficients
    offset = ((1:n_subjects) - (n_subjects + 1) / 2) * .025 * range(xlim);
   
    % plot subject coefficients
    grapeplot(bb+offset,betas(:,bb),...
        'marker','o',...
        'markersize',6.5,...
        'markerfacecolor',subject_clr,...
        'markeredgecolor','w',...
        'linewidth',1.5);
end

% plot animal-pool coefficients
p = stem(1:n_coefficients,bigbetas,...
    'color','k',...
    'marker','o',...
    'markersize',10,...
    'markerfacecolor','k',...
    'markeredgecolor','w',...
    'linewidth',1.5);
p.BaseLine.LineWidth = p.LineWidth;

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end