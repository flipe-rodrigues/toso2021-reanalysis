%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% plot firing rate of ramping & non-ramping neurons across task epochs

% figure initialization
fig = figure(figopt,...
    'name','ramp_firingrate_across_epochs');

% epoch settings
epochs = fieldnames(MeanFR);
n_epochs = numel(epochs);

% axes initialization
xxtick = unique((1:n_epochs)+[-1;0;1]*.05*n_epochs);
xxticklabel = num2cell(xxtick);
xxticklabel(~ismember(xxtick,1:n_epochs)) = {''};
xxticklabel(ismember(xxtick,1:n_epochs)) = strrep(epochs,'_',' ');
axes(axesopt.default,...
    'plotboxaspectratio',[2,1,1],...
    'color','none',...
    'ticklength',axesopt.default.ticklength*.58,...
    'xlim',[1,n_epochs]+[-1,1]*.75,...
    'xtick',xxtick,...
    'xticklabel',xxticklabel,...
    'ylimspec','tight',...
    'clipping','off',...
    'layer','bottom');
ylabel('Firing rate [Hz]');

% preallocation
fr_distro = struct();
fr_avg = struct();
fr_err = struct();

% choice of average and error functions
avgfun = @(x) nanmean(x);
errfun = @(x) [1,1] .* nanstd(x);
avgfun = @(x) nanmedian(x);
errfun = @(x) quantile(x,[.25,.75]) - nanmedian(x);

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    [~,uscore_idcs] = regexp(epoch,'_');
    stim = lower(epoch(1:uscore_idcs(1)-1));
    fr_distro.(epoch) = {...
        MeanFR.(epoch)(ramp_idcs.(stim));...
        MeanFR.(epoch)(nonramp_idcs.(stim))};
    fr_avg.(epoch) = [...
        avgfun(MeanFR.(epoch)(ramp_idcs.(stim)));...
        avgfun(MeanFR.(epoch)(nonramp_idcs.(stim)))];
    fr_err.(epoch) = [...
        errfun(MeanFR.(epoch)(ramp_idcs.(stim)));...
        errfun(MeanFR.(epoch)(nonramp_idcs.(stim)))];
end

% table conversions
fr_avg = struct2table(fr_avg,...
    'rownames',{'ramping','non-ramping'});
fr_err = struct2table(fr_err,...
    'rownames',{'ramping','non-ramping'});

% offset between ramps and non-ramps
xoffsets = [-1,1] * .15;

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
%     plot(ee+xoffsets,rho_avg.(epoch),...
%         'color','k',...
%         'linewidth',1.5);
    for rr = 1 : 2
        errorbar(ee+xoffsets(rr),fr_avg.(epoch)(rr),...
            fr_err.(epoch)(rr,1),fr_err.(epoch)(rr,2),...
            'color','k',...
            'marker','o',...
            'markersize',7.5,...
            'markeredgecolor','k',...
            'markerfacecolor',ramp_clrs(rr,:),...
            'linewidth',1.5);
    end
end

% update axes
yylim = [0,25]; % [min(.5,min(ylim)),1];
yytick = yylim;
yyticklabel = num2cell(yytick);
set(gca,...
    'ylim',yylim + [-1,1] * .1 * range(yylim),...
    'ytick',yytick,...
    'yticklabel',yyticklabel);

% iterate through alignments
for ee = 1 : n_epochs
    xx = [-1,1] * .5 / 3 + ee;
    yy = [1,1] * max(yylim);
    plot(xx,yy,...
        'color','k',...
        'linewidth',1.5);
    epoch = epochs{ee};
    [~,pval] = kstest2(fr_distro.(epoch){1},fr_distro.(epoch){2});
    if pval < .01
        test_str = '**';
    elseif pval < .05
        test_str = '*';
    else
        test_str = 'n.s.';
    end
    text(mean(xx),mean(yy)-.025*range(ylim),test_str,...
        'color','k',...
        'fontsize',16,...
        'horizontalalignment','center',...
        'verticalalignment','bottom');
end

%%
return;
figure;
hold on;
histogram(MeanFR.S1_onset(ramp_idcs.s1));
histogram(MeanFR.S1_onset(nonramp_idcs.s1));

figure;
hold on;
histogram(MeanFR.S1_offset(ramp_idcs.s1));
histogram(MeanFR.S1_offset(nonramp_idcs.s1));

figure;
hold on;
histogram(MeanFR.S2_onset(ramp_idcs.s1));
histogram(MeanFR.S2_onset(nonramp_idcs.s1));

figure;
hold on;
histogram(MeanFR.S2_offset(ramp_idcs.s1));
histogram(MeanFR.S2_offset(nonramp_idcs.s1));