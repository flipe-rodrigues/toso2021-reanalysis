function progressreport( iter, n_iter, operation, bar_len )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    bar_len = 30;
end

progress = floor(iter / n_iter * 100);
if iter > 1
    for jj = 1 : (length(operation) + bar_len + 8)
        fprintf('\b');
    end
end
fprintf('|');
for jj = 1 : bar_len
    if progress/(100/bar_len) >= jj
        fprintf('#');
    else
        fprintf('-');
    end
end
fprintf('| ');
if progress == 100
    fprintf([operation, ' ... done.\n']);
else
    fprintf([operation, ' ... ']);
end

end