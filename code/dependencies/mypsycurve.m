function [psy, obs, theta, bias, slope] = ...
    mypsycurve( intervals, choices, lambda, max_iter )
%MYPSYCURVE Computes the psychometric curve that results from applying a
%logistic regression to the behavior data that is outputted by the
%TAFC-timing task. It models the relationship between the durations of the
%delivered stimuli (intervals) and the forced-choice responses of the
%subject.
%
%   PSY = MYPSYCURVE(DATA) returns the psychometric curve PSY (a struct
%   with three fields: x & y) that best fits our input DATA (a struct
%   outputted by the MYPARSEBHV function).
%
%   PSY = MYPSYCURVE(DATA, LAMBDA) uses a regularization factor LAMBDA (0
%   by default) to fit the data.
%
%   PSY = MYPSYCURVE(DATA, LAMBDA, MAX_ITER) limits the number of
%   iterations (500 by default) performed by FMINUNC when minimizing the
%   cost function computed by COSTFUNCTION.
%
%   [PSY, OBS] = MYPSYCURVE(...)
%
%   [PSY, OBS, THETA] = MYPSYCURVE(...)
%
%   [PSY, OBS, THETA, BIAS] = MYPSYCURVE(...)
%
%   [PSY, OBS, THETA, BIAS, SLOPE] = MYPSYCURVE(...)
%
%   See also LOGREG, FMINUNC, COSTFUNCTION, SIGMOID, MYPARSEBHV.

%% Input Validation
narginchk(2, 4);
if ~isvector(intervals) || ~isvector(choices)
    error('The 1st & 2nd input arguments must be vectors.');
end
if nargin < 3
    lambda = 0;
elseif ~isfloat(lambda) || lambda < 0
    error('The 3rd input argument must be a positive float.');
end
if nargin < 4
    max_iter = 1e6;
elseif ~isscalar(max_iter) || max_iter <= 0 || max_iter ~= floor(max_iter)
    error('The 4th input argument must be a positive integer.');     
end

%% Input Preparation
n_trials = length(intervals);
intervals = reshape(intervals, [n_trials, 1]);
choices = reshape(choices, [n_trials, 1]);
nan_idcs = isnan(intervals) | isnan(choices);
intervals = intervals(~nan_idcs);
choices = choices(~nan_idcs);
temp_stimset = nan(8);
temp_stimset(1:length(unique(intervals))) = unique(intervals);
stimset = [.2, .35, .42, .46, .54, .58, .65, .8];
for ii = 1 : 8
    stim = stimset(ii);
    if ~ismember(stim, temp_stimset)
        stimset(ii) = nan;
    end
end

%% Data Points
obs.x = stimset;
for ii = 1 : length(stimset)
    s_idcs = intervals == stimset(ii);
    obs.y(ii) = sum(choices(s_idcs)) / sum(s_idcs);
    obs.count(ii) = sum(s_idcs);
end

%% Logistic Regression
% theta = logreg(intervals, choices, lambda, max_iter);
options = optimset('gradobj', 'on', 'maxiter', 10e3, 'maxfunevals', 10e3);
theta = glmfit(intervals, choices, 'binomial');
% options = optimset('gradobj', 'on', 'maxiter', max_iter, 'maxfunevals', max_iter);
% initial_theta = [min(obs.y), ...
%                  max(obs.y) - min(obs.y), ...
%                  nanmean(stimset), ...
%                  max(obs.y)];
% theta = fminsearch(@(t)sig4sse(t, stimset, obs.y),initial_theta, options);

%% Psychometric Fit
psy.x = linspace(min(stimset),max(stimset),1e4);
psy.y = sigmoid(theta(1) + theta(2) * psy.x);
% psy.y = theta(1) + theta(2) ./ (1 + exp(-1 * theta(3) * (psy.x - theta(4))));

[~, idx] = min(abs(psy.y - .5));
bias = psy.x(idx);
try
    slope = diff(psy.y(idx - 1 : idx + 1)) / diff(psy.x(idx - 1 : idx + 1));
end

%%
function [ error ] = sig4sse(betas, X, Y )
    Yhat = betas(1) + betas(2)./(1+exp(-1*betas(3)*(X-betas(4))));
    error = nansum((Y-Yhat).^2);
end

end