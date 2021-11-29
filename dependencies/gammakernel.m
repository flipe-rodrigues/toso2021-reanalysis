function k = gammakernel(varargin)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
     
    p = inputParser;
    p.addParameter('peakx',10);
    p.addParameter('binwidth',2);
    p.parse(varargin{:});   
    
    k = p.Results;
    k.type = 'gamma';
    k.shape = 2;
    k.scale = k.peakx / k.binwidth;
    k.nbins = round(k.peakx / k.binwidth * 20);
    k.paddx = [-1,1] / 2 * k.nbins * k.binwidth;
    k.bins = -k.nbins / 2 + 1 : k.nbins / 2;
    k.x = k.bins * k.binwidth;
    k.pdf = gampdf(k.bins,k.shape,k.scale);
    k.pdf = k.pdf / sum(k.pdf);
end

