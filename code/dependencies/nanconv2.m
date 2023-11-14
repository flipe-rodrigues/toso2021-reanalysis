function C = nanconv2(A,x,y)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    
    % input parsing
    [n_rows,n_cols] = size(A);
    n = length(x);
    m = length(y);
    
    %
    x = flipud(x(:));
    y = flipud(y(:));

    % preallocation
    C = double(A);
    
    %
    B1 = padarray(C,floor([n,1]/2),nan,'both');

    %
    for ii = 1 : n_rows
        row_idcs = (1 : n) + ii - 1;
        nan_flags = all(isnan(B1(row_idcs,:)),2);
        if all(nan_flags) || nansum(x(~nan_flags)) == 0
            continue;
        end
        u = x(~nan_flags) / nansum(x(~nan_flags));
        for jj = 1 : n_cols
            b = B1(row_idcs,jj);
            C(ii,jj) = u' * b(~nan_flags);
        end
    end
    
    %
    B2 = padarray(C,floor([1,m]/2),nan,'both');
    
    %
    for jj = 1 : n_cols
        col_idcs = (1 : m) + jj - 1;
        nan_flags = all(isnan(B2(:,col_idcs)),1);
        if all(nan_flags) || nansum(y(~nan_flags)) == 0
            continue;
        end
        v = y(~nan_flags) / nansum(y(~nan_flags));
        for ii = 1 : n_rows
            b = B2(ii,col_idcs)';
            C(ii,jj) = v' * b(~nan_flags);
        end
    end
end

