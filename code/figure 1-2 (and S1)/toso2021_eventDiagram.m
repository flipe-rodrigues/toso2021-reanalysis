%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% seed
rng(0);

%% vibration generation settings
fs = 10e3;              % sampling frequency
dt = 1 / fs;            % sampling period
downsampling_factor = 50;

%% trace settings

% time settings
m = 1e3;
gt_eg = 100;
rt_eg = 300;
mt_eg = 600;
pre_t1_eg = unique(pre_s1_delay(valid_flags));
t1_eg = t1_set(t1_mode_idx);
i1_eg = i1_set(i1_mode_idx);
t2_eg = max(unique(t2(valid_flags & t1 == t1_eg)));
i2_eg = max(unique(i2(valid_flags & i1 == i1_eg)));
post_t2_eg = post_s2_delay;
trial_dur = pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg;
lt_eg = trial_dur * .165;
t = linspace(-lt_eg,trial_dur + rt_eg + mt_eg + 200,m);

% sample stimulus speeds
s1_m = round(t1_eg / 1e3 * fs / downsampling_factor);
s1_t = linspace(0,t1_eg,s1_m) + pre_t1_eg;
s1_eg = abs(normrnd(0,i1_eg,s1_m,1));
s1_eg = (s1_eg - min(s1_eg)) / range(s1_eg);
s1_eg([1,end]) = 0;
s2_m = round(t2_eg / 1e3 * fs / downsampling_factor);
s2_t = linspace(0,t2_eg,s2_m) + pre_t1_eg + t1_eg + isi;
s2_eg = abs(normrnd(0,i2_eg,s2_m,1));
s2_eg = (s2_eg - min(s2_eg)) / range(s2_eg);
s2_eg([1,end]) = 0;

% sample go cue
go_m = 1e3;
go_t = linspace(0,gt_eg,go_m) + pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg + 3;
go_f = 2.5 * 1e3 / gt_eg;
go_eg = sin(2 * pi * linspace(0,gt_eg,go_m) / 1e3 * go_f) * 1/2;
go_eg([1,end]) = 0;

% preallocation
traces = struct();
traces.init = zeros(m,1);
traces.s1 = zeros(m,1);
traces.s2 = zeros(m,1);
traces.go = zeros(m,1);
traces.choice = zeros(m,1);
traces.rwd = zeros(m,1);
traces = struct2table(traces);
n_traces = size(traces.Variables,2);
trace_lbls = traces.Properties.VariableNames;

% text annotations
trace_str.init = {'Initiation port'};
trace_str.s1 = {'Stimulus 1 (S1)'};
trace_str.s2 = {'Stimulus 2 (S2)'};
trace_str.go = {'Go cue'};
trace_str.choice = {'Choice spout (T2 > T1)'};
trace_str.rwd = {'Reward'};

% initiation port trace
init_flags = t >= 0 & ...
    t <= trial_dur + rt_eg;
traces.init(init_flags) = 1;

% s1 trace
s1_flags = t >= pre_t1_eg & ...
    t <= pre_t1_eg + t1_eg;
traces.s1(s1_flags) = nan;

% s2 trace
s2_flags = t >= pre_t1_eg + t1_eg + isi & ...
    t <= pre_t1_eg + t1_eg + isi + t2_eg;
traces.s2(s2_flags) = nan;

% go trace
go_flags = t >= pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg & ...
    t <= pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg + gt_eg;
traces.go(go_flags) = 1;

% choice spout trace
choice_flags = t >= pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg + rt_eg + mt_eg;
traces.choice(choice_flags) = 1;

% reward trace
rwd_flags = t >= pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg + rt_eg + mt_eg & ...
    t <= pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg + rt_eg + mt_eg + gt_eg;
traces.rwd(rwd_flags) = 1;

%% plot event diagram

% figure initialization
fig = figure(figopt,...
    'position',[279,240,980,520],...
    'name','event_diagram',...
    'color','w');

% axes initialization
aspect = [2,1,1];
set(gca,axesopt.default,...
    'xlim',[min(t),max(t)],...
    'ylim',[-n_traces,0]+[-1,1]*n_traces*.05*aspect(1)/aspect(2),...
    'plotboxaspectratio',aspect,...
    'color','none',...
    'xcolor','none',...
    'ycolor','none',...
    'clipping','off',...
    'layer','top');

% plot settings
linewidth = 1;

% iterate through traces
for tt = 1 : n_traces
    trace_lbl = trace_lbls{tt};
    offset = -tt;
    gain = 1/3;
    stairs(t,traces.(trace_lbl) * gain + offset,...
        'color','k',...
        'linewidth',linewidth);
    
    % plot signals
    if ismember(trace_lbl,{'s1','s2'})
        plot(eval([trace_lbl,'_t']),eval([trace_lbl,'_eg']) * 1/2 + offset,...
            'color','k',...
            'linewidth',linewidth);
    end
    
    % text annotations
    text(min(xlim),offset+.05,trace_str.(trace_lbl),...
        'color','k',...
        'fontsize',12,...
        'horizontalalignment','left',...
        'verticalalignment','bottom');
end

% reference lines
plot([1,1]*0,[-n_traces-1*gain,0],':k');
plot([1,1]*pre_t1_eg,[-n_traces-1*gain,0],':k');
plot([1,1]*pre_t1_eg+t1_eg,[-n_traces-1*gain,0],':k');
plot([1,1]*pre_t1_eg+t1_eg+isi,[-n_traces-1*gain,0],':k');
plot([1,1]*pre_t1_eg+t1_eg+isi+t2_eg,[-n_traces-1*gain,0],':k');
plot([1,1]*pre_t1_eg+t1_eg+isi+t2_eg+post_t2_eg,[-n_traces-1*gain,0],':k');
plot([1,1]*pre_t1_eg+t1_eg+isi+t2_eg+post_t2_eg+rt_eg,[-n_traces-1*gain,0],':k');
plot([1,1]*pre_t1_eg+t1_eg+isi+t2_eg+post_t2_eg+rt_eg+mt_eg,[-n_traces-1*gain,0],':k');

% delay lines
delay_h = struct();
delay_h.pres1 = plot([0,pre_t1_eg],[1,1]*-1.15,...
    'color',[1,1,1]*.65,...
    'linewidth',linewidth);
delay_h.s1 = plot(pre_t1_eg+[0,t1_eg],[1,1]*-2.15,...
    'color',[1,1,1]*.0,...
    'linewidth',linewidth);
delay_h.isi = plot(pre_t1_eg+t1_eg+[0,isi],[1,1]*-2.15,...
    'color',[1,1,1]*.65,...
    'linewidth',linewidth);
delay_h.s2 = plot(pre_t1_eg+t1_eg+isi+[0,t2_eg],[1,1]*-3.15,...
    'color',[1,1,1]*.0,...
    'linewidth',linewidth);
delay_h.posts2 = plot(pre_t1_eg+t1_eg+isi+t2_eg+[0,post_t2_eg],[1,1]*-3.15,...
    'color',[1,1,1]*.65,...
    'linewidth',linewidth);
delay_h.rt = plot(pre_t1_eg+t1_eg+isi+t2_eg+post_t2_eg+[0,rt_eg],[1,1]*-4.15,...
    'color',[1,1,1]*.65,...
    'linewidth',linewidth);
delay_h.mt = plot(pre_t1_eg+t1_eg+isi+t2_eg+post_t2_eg+rt_eg+[0,mt_eg],[1,1]*-5.15,...
    'color',[1,1,1]*.65,...
    'linewidth',linewidth);

% delay strings
delay_s = struct();
delay_s.pres1 = 'Pre-S1 delay';
delay_s.s1 = 'T1';
delay_s.isi = 'Inter-stimulus interval (ISI)';
delay_s.s2 = 'T2';
delay_s.posts2 = 'Post-S2 delay';
delay_s.rt = 'Reaction time';
delay_s.mt = 'Movement time';

% delay labels
delay_fields = fieldnames(delay_h);
n_delays = numel(delay_fields);
for dd = 1 : n_delays
    xpos = mean(delay_h.(delay_fields{dd}).XData);
    ypos = mean(delay_h.(delay_fields{dd}).YData);
    str = delay_s.(delay_fields{dd});
    text(xpos,ypos - .005 * range(ylim),str,...
        'color',delay_h.(delay_fields{dd}).Color,...
        'fontsize',10,...
        'horizontalalignment','center',...
        'verticalalignment','top');
end

text(pre_t1_eg+t1_eg+100,-2 + .05,'I1',...
    'color',delay_h.s1.Color,...
    'fontsize',10,...
    'horizontalalignment','center',...
    'verticalalignment','bottom');
text(pre_t1_eg+t1_eg+isi+t2_eg+100,-3 + .05,'I2',...
    'color',delay_h.s2.Color,...
    'fontsize',10,...
    'horizontalalignment','center',...
    'verticalalignment','bottom');

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end