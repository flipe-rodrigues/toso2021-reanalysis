%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% parse session data

% flag (unique) session transitions
session_start_idcs = find(data.Trial == 1);
valid_session_flags = ...
    ismember(session_start_idcs,find(valid_flags & unique_flags));
session_start_idcs = session_start_idcs(valid_session_flags);
n_total_sessions = numel(session_start_idcs);
session_bound_idcs = unique([session_start_idcs;n_total_trials]);

% preallocation
session_idcs = nan(n_total_trials,1);
session_trial_count = nan(n_total_sessions,1);
session_neuron_count = nan(n_total_sessions,1);
session_subject = nan(n_total_sessions,1);

% iterate through sessions
for ss = 1 : n_total_sessions
    progressreport(ss,n_total_sessions,'parsing session data');
    idcs = session_bound_idcs(ss) : session_bound_idcs(ss+1) - 1;
    session_idcs(idcs) = ss;
    session_trial_count(ss) = sum(...
        valid_flags & ...
        unique_flags & ...
        ismember(1:n_total_trials,idcs)');
    session_neurons = unique(data.NeuronNumb(...
        valid_flags & ...
        ismember(1:n_total_trials,idcs)' & ...
        ismember(data.NeuronNumb,flagged_neurons)));
    session_neuron_count(ss) = numel(session_neurons);
    session_subject(ss) = unique(subjects(...
        valid_flags & ...
        unique_flags & ...
        ismember(1:n_total_trials,idcs)'));
end