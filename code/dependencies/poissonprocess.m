function [ count, times ] = poissonprocess( lambda, T )
%POISSONPROCESS Draws events from a poisson process.
%   [COUNT] = POISSONPROCESS(LAMBDA,T) returns the number of events drawn
%   from a poisson process parameterized by the success rate LAMBDA over
%   an interval of duration T. If LAMBDA is a scalar, then a constant rate
%   is assumed (meaning the poisson process is homogeneous). If LAMBDA is a
%   vector, then an inhomogenous process is assumed, wherein the success
%   rate changes over time T as defined by LAMBDA.
%
%   [COUNT,TIMES] = POISSONPROCESS(LAMBDA,T) returns the time of each
%   occurrence in the column vector TIMES, in addition to the event COUNT.
%
%   Examples:
%
%      Example 1: Calculate the number of events that happen over a period
%      of 5 seconds in a homogeneous poisson process with a rate of 10 Hz.
%         n = poissonprocess(10,5);
%
%      Example 2: Draw event times from a homogenous poisson process with
%      a success rate of 10 Hz and a duration of 5 seconds.
%         [~,t] = poissonprocess(10,5);
%
%      Example 3: Calculate the number of events and their respective time
%      of occurrence from an inhomogeneous poisson process with a
%      time-varying rate that increases from 1 to 10 Hz over 5 seconds.
%         [n,t] = poissonprocess(1:10,5);
%
%
%   See also POISSRND, UNIFRND.

%   Author: Filipe S. Rodrigues.

narginchk(2,2);

m = length(lambda);

% homogenous poisson process
if isscalar(lambda)
    count = poissrnd(lambda * T);
    times = sort(unifrnd(0,T,count,1));

% inhomogenous poisson process
else
    maxlambda = max(lambda);
    count = poissrnd(maxlambda * T);
    times = sort(unifrnd(0,T,count,1));
    for ii = 1 : count
        die = unifrnd(0,1);
        idx = ceil(times(ii) / T * m);
        if die > (lambda(idx) / maxlambda)
            count = count - 1;
            times(ii) = nan;
        end
    end
    times = times(~isnan(times));
end

end