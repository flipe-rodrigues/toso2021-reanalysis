function p = errorpatch(x,y,e,facecolor,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
p.addParameter('facealpha',1);
p.addParameter('edgecolor','none');
p.addParameter('parent',gca);
p.parse(varargin{:});
patchopt = p.Results;

if iscolumn(x)
    x = x';
end
if iscolumn(y)
    y = y';
end
if iscolumn(e)
    e = e';
end

nan_flags = isnan(x) | isnan(y) | isnan(e);

X = [x(~nan_flags),fliplr(x(~nan_flags))];
Y = [y(~nan_flags)-e(~nan_flags),fliplr(y(~nan_flags)+e(~nan_flags))];

patch(X,Y,facecolor,patchopt);

end

