%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% plot stereotypy of ramping & non-ramping neurons across task epochs

% figure initialization
fig = figure(figopt,...
    'name','ramp_proportion_across_epochs');

% epoch settings
epochs = fieldnames(ramp_idcs_ud);
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
    'ylim',[0,1]+[-1,1]*.0,...
    'ytick',linspace(0,1,5),...
    'clipping','off',...
    'layer','bottom');
ylabel('Proportion of neurons');

% preallocation
proportion_ramp = struct();
proportion_stereo = struct();
proportion_stim = struct();

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    proportion_ramp.(epoch) = [...
        numel(ramp_idcs_ud.(epoch)); numel(nonramp_idcs_ud.(epoch))] ./ ...
        (numel(ramp_idcs_ud.(epoch)) + numel(nonramp_idcs_ud.(epoch)));
    proportion_stereo.(epoch) = [...
        numel(stereo_idcs_ud.(epoch)); numel(nonstereo_idcs_ud.(epoch))] ./ ...
        (numel(stereo_idcs_ud.(epoch)) + numel(nonstereo_idcs_ud.(epoch)));
end

%
proportion_stim.s1 = [...
    numel(ramp_idcs.s1); numel(nonramp_idcs.s1)] ./ ...
    (numel(ramp_idcs.s1) + numel(nonramp_idcs.s1));
proportion_stim.s2 = [...
    numel(ramp_idcs.s2); numel(nonramp_idcs.s2)] ./ ...
    (numel(ramp_idcs.s2) + numel(nonramp_idcs.s2));

% table conversions
proportion_ramp = struct2table(proportion_ramp,...
    'rownames',{'ramping','non-ramping'});
proportion_stereo = struct2table(proportion_stereo,...
    'rownames',{'stereotypycal','non-stereotypycal'});
proportion_stim = struct2table(proportion_stim,...
    'rownames',{'ramping','non-ramping'});

% offset between ramps and non-ramps
xoffsets = [-1,1] * .15;

% iterate through alignments
for ee = 1 : n_epochs
    epoch = epochs{ee};
    xxrange = xxtick(ee*3+[-2:0]);
    for rr = 1 : 2
        xx = xxrange(rr*2-1) + [0,.05] * n_epochs * (-1)^(rr == 2) * .85;
        xpatch = [xx,fliplr(xx)];
        ypatch = sort(repmat([0,proportion_ramp.(epoch)(rr)],1,2));
        patch(xpatch,ypatch,ramp_clrs(rr,:),...
            'facealpha',1,...
            'edgecolor','k',...
            'linewidth',1.5);
    end
end

%%
return;
figure;
hold on;
histogram(StereoCrit.S1_onset(ramp_idcs.s1));
histogram(StereoCrit.S1_onset(nonramp_idcs.s1));

figure;
hold on;
histogram(StereoCrit.S1_offset(ramp_idcs.s1));
histogram(StereoCrit.S1_offset(nonramp_idcs.s1));

figure;
hold on;
histogram(StereoCrit.S2_onset(ramp_idcs.s1));
histogram(StereoCrit.S2_onset(nonramp_idcs.s1));

figure;
hold on;
histogram(StereoCrit.S2_offset(ramp_idcs.s1));
histogram(StereoCrit.S2_offset(nonramp_idcs.s1));