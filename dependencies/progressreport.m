function progressreport( iter, n_iter, message, bar_len )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 4
        bar_len = 30;
    end
    persistent progress;
    if isempty(progress)
        progress = -1;
    elseif progress ~= floor(iter / n_iter * 100)
        progress = floor(iter / n_iter * 100);
    else
       return; 
    end
    if iter > 1
        for jj = 1 : (length(message) + bar_len + 8)
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
        fprintf([message, ' ... done.\n']);
    else
        fprintf([message, ' ... ']);
    end
end