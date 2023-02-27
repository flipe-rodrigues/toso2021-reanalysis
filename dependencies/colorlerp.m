function map = colorlerp(clrs,k)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    K = 1e3;
    [n,d] = size(clrs);
    m = K / (n - 1);
    map = nan(k,d);
    for ii = 1 : n - 1
        idcs = (1 : m) + (ii - 1) * m;
        map(idcs,:) = cell2mat(...
            arrayfun(@(dim)linspace(clrs(ii,d),clrs(ii+1,d),m)',1:d,...
            'uniformoutput',false));
    end
    idcs2keep = round(linspace(1,K,k));
    idcs2keep = min(max(idcs2keep,1),K);
    map = map(idcs2keep,:);
end

