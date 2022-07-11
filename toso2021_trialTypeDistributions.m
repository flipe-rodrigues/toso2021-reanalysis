%% initialization
if ~exist('data','var')
    toso2021_wrapper;
end

%% construct S2-aligned, Ti- & Ii-split psths

% time settings
roi = [0,t_set(end)];
roi_n_bins = range(roi) * psthbin;
roi_time = linspace(roi(1),roi(2),roi_n_bins);

% preallocation
trial_type_numbers = nan(numel(neuron_idcs),n_stimuli,n_contrasts);

% iterate through neurons
for nn = neuron_idcs'
    progressreport(nn,numel(neuron_idcs),'parsing neural data');
    neuron_flags = data.NeuronNumb == neuron_idcs(nn);
    
    % iterate through stimuli
    for tt = 1 : n_stimuli
        stim_flags = stimuli == stim_set(tt);
        
        % iterate through contrasts
        for ii = 1 : n_contrasts
            contrast_flags = contrasts == contrast_set(ii);
            s2_spike_flags = ...
                valid_flags & ...
                neuron_flags & ...
                stim_flags & ...
                contrast_flags;
            if sum(s2_spike_flags) == 0
                continue;
            end
            
            % fetch spike counts & compute spike rates
            s2_spike_counts = data.FR(s2_spike_flags,:);
            s2_spike_rates = conv2(...
                1,kernel.pdf,s2_spike_counts,'valid')' / psthbin * 1e3;
            s2_n_trials = size(s2_spike_counts,1);
            
            % T2-aligned spike rates
            s2_alignment_offset = ...
                pre_init_padding + ...
                pre_t1_delay(s2_spike_flags) + ...
                t1(s2_spike_flags) + ...
                isi;
            s2_alignment_flags = ...
                valid_time >= s2_alignment_offset + roi(1) & ...
                valid_time < s2_alignment_offset + t2(s2_spike_flags);
            s2_chunk_flags = ...
                valid_time >= s2_alignment_offset + roi(1) & ...
                valid_time < s2_alignment_offset + roi(2);
            s2_spkrates = s2_spike_rates;
            s2_spkrates(~s2_alignment_flags') = nan;
            s2_spkrates = reshape(...
                s2_spkrates(s2_chunk_flags'),[roi_n_bins,s2_n_trials])';
            
            % neuron selection criteria
            trial_type_numbers(nn,tt,ii) = s2_n_trials;
        end
    end
end

%% plot trial type distributions

% figure & axes initialization
fig = figure(figopt,...
    'name',sprintf('trial_type_distributions_%s',contrast_str));
axes(axesopt.default,...
    'ylim',[1,n_stimuli] + [-1,1] * .05 * n_stimuli,...
    'ytick',1:n_stimuli,...
    'yticklabel',num2cell(stim_set),...
    'xscale','log');
xlabel('Trial count');
ylabel(stim_lbl);
zlabel('Neuron count');
view(30,80);

% bin settings
n_bins = 50;
x_edges = linspace(0,50,n_bins+1);
y_edges = linspace(0,n_stimuli+1,n_bins+1);

% iterate through contrasts
for kk = 1 : n_contrasts

    % iterate through stimuli
    for ii = 1 : n_stimuli
        
        % compute trial counts
        bin_counts = histcounts2(squeeze(trial_type_numbers(:,ii,kk)),...
            repmat(ii+kk/10,size(trial_type_numbers,1),1),...
            'xbinedges',x_edges,...
            'ybinedges',y_edges);
        
        % plot trial counts
        histogram2('xbinedges',x_edges,...
            'ybinedges',y_edges,...
            'bincounts',bin_counts,...
            'edgecolor','none',...
            'facecolor',contrast_clrs(kk,:));
    end
end

% save figure
if want2save
    svg_file = fullfile(panel_path,[fig.Name,'.svg']);
    print(fig,svg_file,'-dsvg','-painters');
end