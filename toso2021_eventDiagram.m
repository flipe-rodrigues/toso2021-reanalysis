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
traces = struct2table(traces);
n_traces = size(traces.Variables,2);
trace_lbls = traces.Properties.VariableNames;

% text annotations
trace_str.init = {'Initiation port'};
trace_str.s1 = {'Stimulus 1 (S1)'};
trace_str.s2 = {'Stimulus 2 (S2)'};
trace_str.go = {'Go cue'};
trace_str.choice = {'Choice spout'};

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
traces.go(go_flags) = nan;

% choice spout trace
choice_flags = t >= pre_t1_eg + t1_eg + isi + t2_eg + post_t2_eg + rt_eg + mt_eg;
traces.choice(choice_flags) = 1;

%% plot event diagram

% figure initialization
fig = figure(figopt,...
    'position',[15,100,1000,520],...
    'name','event_diagram',...
    'color','w');

% axes initialization
aspect = [2.25,1,1];
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
    if ismember(trace_lbl,{'s1','s2','go'})
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

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end

return;

%% vibration generation settings
fs = 10e3;              % sampling frequency
dt = 1 / fs;            % sampling period
m = 100e3;              % length of signal
t = (0 : m - 1) * dt;   % time vector

% butterworth filter
fc = 150;
[b,a] = butter(2,fc/(fs/2));

%% trace settings

%
mus = zeros(n_i,1);
sigs = i_set ./ sqrt(2 / pi);

% preallocation
v_samples = nan(n_i,m);

figure; hold on;
title('Whisker vibration velocity distributions')
xlabel('Velocity (mm.s^{-1})')
ylabel('PDF')

% iterate through intensities
for ii = n_i : -1 : 1
    v = normrnd(mus(ii),sigs(ii),1,m); % + sin(2 * pi * t * i_set(ii)) * 10;
    v_samples(ii,:) = v; % filter(b,a,v);
    [mean(abs(v)),mean(abs(v_samples(ii,:))),i_set(ii)]
    histogram(v_samples(ii,:),linspace(-750,750,100),...
        'facecolor',i2_clrs(ii,:),...
        'normalization','pdf');
end

figure; hold on;
title('Whisker vibration speed distributions')
xlabel('Speed (mm.s^{-1})')
ylabel('PDF')

% iterate through intensities
for ii = n_i : -1 : 1
    histogram(abs(v_samples(ii,:)),linspace(0,750,100),...
        'facecolor',i2_clrs(ii,:),...
        'normalization','pdf');
end

%%
figure; hold on;
title('Single-Sided Amplitude Spectrum of v(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
% set(gca,'yscale','log');
xlim([0,500]);

% iterate through intensities
for ii = n_i : -1 : 1
    X = v_samples(ii,:);
    Y = fft(X);
    P2 = abs(Y/m);
    P1 = P2(1:m/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(m/2))/m;
    plot(f,P1,...
        'color',i2_clrs(ii,:))
end

%%

figure;

% iterate through intensities
for ii = n_i : -1 : 1
    subplot(1,n_i,ii);
    spectrogram(v_samples(ii,:),256,250,256,fs,'yaxis');
    axis tight
end