function img = tensor2rgb(tensor, clrs, brightness)

if nargin < 3
    brightness = 1;
end

% input parsing
[width,height,dim] = size(tensor);

% normalization
tensor = (tensor - min(tensor,[],'all')) ./ range(tensor,'all');

% transformation matrix
A = (1 - clrs) / eye(3);

% color conversion
pixels = brightness - reshape(tensor,width*height,dim) * A;

% reshape from pixels to image
img = reshape(pixels,width,height,3);

% normalization
% img = (img - min(img,[],'all')) ./ range(img,'all');

end

