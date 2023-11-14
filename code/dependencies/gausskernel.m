function k = gausskernel(varargin)
    %GAUSSKERNEL Creates a gaussian-shaped smoothing kernel
    %   K = GAUSSKERNEL('sig',SIG,'binwidth',BINWIDTH) creates a gaussian
    %   kernel with a standard deviation of SIG and a resolution of BINWIDTH.
     
    p = inputParser;
    p.addParameter('sig',10);
    p.addParameter('binwidth',2e-3);
    p.parse(varargin{:});   
    
    k = p.Results;
    k.nbins = 2 * round(k.sig / k.binwidth * 9 / 2);
    k.paddx = [-1,1] / 2 * k.nbins * k.binwidth;
    k.bins = -k.nbins / 2 + 1 : k.nbins / 2;
    k.nbins = numel(k.bins);
    k.x = k.bins * k.binwidth;
    k.pdf = normpdf(k.x,0,k.sig);
    k.pdf = k.pdf / sum(k.pdf);
end

