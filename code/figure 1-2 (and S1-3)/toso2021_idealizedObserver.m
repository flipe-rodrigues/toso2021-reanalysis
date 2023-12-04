%% check 'main.m' has run (and run it if not)
toso2021_maincheck;

%% idealized observer settings
beta_nsd = 3;
beta_ndd = 1;
% s1 = datasample(s_set,n_total_trials);
% s2 = datasample(s_set,n_total_trials);
% d1 = datasample(d_set,n_total_trials);
% d2 = datasample(d_set,n_total_trials);
% nsd = round((s2 - s1) ./ (s2 + s1),2);
% ndd = round((d2 - d1) ./ (d2 + d1),1);
p_choice = 1 ./ (1 + exp(-beta_nsd * nsd - beta_ndd * ndd));

%
% tfun = @(x) log(x);
% invfun = @(x) exp(x);
% si_bounds = tfun([s_set(1),s_set(end)]) + [-1,1] * .1 * range(tfun(s_set));
% di_bounds = tfun([d_set(1),d_set(end)]) + [-1,1] * .1 * range(tfun(d_set));
% si_x = linspace(invfun(si_bounds(1)),invfun(si_bounds(2)),1e2);
% di_x = linspace(invfun(di_bounds(1)),invfun(di_bounds(2)),1e2);
% p_choice_mat = 1 ./ (1 + exp(...
%     -beta_nsd * (si_x' - si_x) ./ (si_x' + si_x) + ...
%     -beta_ndd * (di_x' - di_x) ./ (di_x' + di_x)));

%% simulated idealized observer
if ~exist('idealized_choice','var')
    idealized_choice = rand(n_total_trials,1) < p_choice;
end