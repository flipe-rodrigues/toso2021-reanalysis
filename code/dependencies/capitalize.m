function [str] = capitalize(str)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    str = sprintf('%s%s',upper(str(1)),lower(str(2:end)));
end

