function [s] = sigmoid( z )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

s = 1 ./ (1 + exp(-z));

end

