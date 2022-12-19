%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% trace settings

% 

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