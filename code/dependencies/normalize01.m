function z = normalize01(x,d)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    d = 1 + isrow(x);
end
z = (x - min(x,[],d)) ./ range(x,d);
end

