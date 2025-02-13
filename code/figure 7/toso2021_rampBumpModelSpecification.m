%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% ROI settings
roi = [0,t_set(end)];
T = range(roi) / psthbin;
ti = roi(1);
tf = roi(2);
t = linspace(ti,tf,T);
t_units = 1e3;
dt = diff(t(1:2));

%% model specification

% figure initialization
fig = figure(figopt,...
    'position',[250,630,415,420],...
    'name','ramps_model_specification',...
    'color',[1,1,1]*245/255);

% model spec
mu_mdl = .5 * range(roi);
gamma_mdl = 5;
sigma_mdl = .15 * range(roi);
x_mdl = generativerate(t,gamma_mdl,mu_mdl,0,sigma_mdl);

% axes initialization
axes(axesopt.default,...
    'xcolor','k',...
    'xlim',roi+[-1,1]*.05*range(roi),...
    'xtick',unique([roi,mu_mdl+[-1,0,1]*sigma_mdl]),...
    'xticklabel',{num2str(roi(1)),'$\mu-\sigma$','$\mu$','$\mu+\sigma$',num2str(roi(2))},...
    'ylim',[0,gamma_mdl]+[-1,1]*.15*range([0,gamma_mdl]),...
    'ytick',[0,gamma_mdl],...
    'yticklabel',{'$b$','$\gamma$'},...
    'ticklabelinterpreter','latex',...
    'clipping','off',...
    'plotboxaspectratio',[3,1,1]);
xlabel('Time (ms)');
ylabel('Firing rate (Hz)');

% model model illustration
plot(t,x_mdl,...
    'color','k',...
    'linewidth',1.5);

% annotate model spcification
text(.5,1.1,...
    ['$\lambda_{n,k}(t)\propto\gamma_{n}\cdot\mathcal{N}',...
    '(\nu_{k}\sim\mathcal{N}(\mu_{n},\eta_{n}),\sigma_{n})$'],...
    'color','k',...
    'fontsize',12,...
    'interpreter','latex',...
    'horizontalalignment','center',...
    'units','normalized');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

%% example neurons with opposite gammas & etas
K_eg = 10;

% preallocation
N_eg = 2;
X_eg = nan(T,N_eg,K_eg);
R_eg = nan(T,N_eg,K_eg);

% parameter choices
mus_eg = [...
    [1,1].*.75;...
    [1,1].*.25] * tf;
etas_eg = [...
    [0,0];...
    [0,1]*.25] * tf;
gammas_eg = [...
    [4,2];...
    [3,3]];
sigma_eg = .15 * (tf - ti);
n_eg = size(mus_eg,1);

% iterate through examples
for ee = 1 : n_eg
    
    % figure initialization
    fig = figure(figopt,...
        'position',[250+(ee-1)*215,350,200,275],...
        'name',sprintf('ramps_simulation_eg%i',ee),...
        'color',[1,1,1]*245/255);
    
    % axes settings
    yspacer = mean(gammas_eg(ee,:));
    axes(axesopt.default,...
        'xlim',[ti,tf],...
        'xtick',[ti,tf],...
        'ytick',[1,K_eg*yspacer],...
        'ytick',(1:K_eg)*yspacer,...
        'yticklabel',num2cell(1:K_eg),...
        'ylim',[1,K_eg+1]*yspacer+[-1,0]*.05*K_eg*yspacer,...
        'ylimspec','tight',...
        'ycolor','none',...
        'clipping','off',...
        'ticklength',axesopt.default.ticklength*2,...
        'plotboxaspectratio',[1,1,1]);
    xlabel('Time (ms)');
    ylabel('Trials');
    
    % iterate through neurons
    for nn = 1 : N_eg
        
        % iterate through trials
        for kk = 1 : K_eg
            X_eg(:,nn,kk) = generativerate(...
                t,gammas_eg(ee,nn),mus_eg(ee,nn),etas_eg(ee,nn),sigma_eg);
            [n,ts] = poissonprocess(X_eg(:,nn,kk),(tf-ti)/t_units);
            spk_times = ts * t_units + ti;
            spk_counts = histcounts(spk_times,...
                'binedges',[t,t(end)+1]);
            R_eg(:,nn,kk) = movsum(spk_counts,dt) / (dt/t_units);
        end
    end
    
    % compute cross-trial mean
    x_eg = nanmean(X_eg,3);
    
    % iterate through trials
    for kk = K_eg : -1 : 1
        
        % iterate through neurons
        for nn = 1 : N_eg
            
            % plot single trial
            plot(t,X_eg(:,nn,kk)+kk*yspacer,...
                'color','w',...
                'linewidth',1.5);
            xpatch = [fliplr(t),t];
            ypatch = [zeros(T,1);X_eg(:,nn,kk)] + kk * yspacer;
            patch(xpatch,ypatch,ramp_clrs(nn,:),...
                'facealpha',1,...
                'edgecolor','none',...
                'linewidth',1);
        end
    end
    
    % parameter annotation
    if ee == 1
        text(.5,1.3,'$\mu_1=\mu_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.2,'$\gamma_1>\gamma_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.1,'$\eta_1=\eta_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
    else
        text(.5,1.3,'$\mu_1=\mu_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.2,'$\gamma_1=\gamma_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
        text(.5,1.1,'$\eta_1>\eta_2$',...
            'color','k',...
            'interpreter','latex',...
            'fontsize',12,...
            'horizontalalignment','center',...
            'units','normalized');
    end
    
    % save figure
    if want2save
        svg_file = fullfile(panel_path,[fig.Name,'.svg']);
        print(fig,svg_file,'-dsvg','-painters');
    end
end