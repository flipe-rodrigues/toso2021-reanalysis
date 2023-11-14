function map = colorlerp(clrs,k)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    K = 1e3;
    n = size(clrs,1);
    m = K / (n - 1);
    map = nan(k,3);
    for ii = 1 : n - 1
        idcs = (1 : m) + (ii - 1) * m;
        map(idcs,:) = [...
            linspace(clrs(ii,1),clrs(ii+1,1),m); ...
            linspace(clrs(ii,2),clrs(ii+1,2),m); ...
            linspace(clrs(ii,3),clrs(ii+1,3),m); ...
            ]';
    end
    idcs2keep = round(linspace(1,K,k));
    idcs2keep = min(max(idcs2keep,1),K);
    map = map(idcs2keep,:);
end

