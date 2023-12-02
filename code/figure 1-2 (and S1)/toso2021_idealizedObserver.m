%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% idealized observer settings
beta_nsd = 2.5;
beta_ndd = 1;
p_choice = 1 ./ (1 + exp(-beta_nsd * nsd - beta_ndd * ndd));

%% simulated idealized observer
if ~exist('idealized_choice','var')
    idealized_choice = rand(n_total_trials,1) < p_choice;
end