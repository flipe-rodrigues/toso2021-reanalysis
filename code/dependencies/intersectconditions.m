function conditions = intersectconditions(varargin)
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    % output preallocation
    conditions = struct();

    % input parsing
    feature_names = varargin(1:4:end);
    feature_values = [varargin{2:4:end}];
    inclusion_constraints = varargin(3:4:end);
    exclusion_contraints = varargin(4:4:end);
    [n_trials,n_features] = size(feature_values);
    feature_flags = false(n_trials,n_features);

    %
    feature_sets = cell(n_features,1);
    for ff = 1 : n_features
        complete_set = unique(feature_values(:,ff))';
        feature_sets{ff} = ...
            complete_set(~ismember(complete_set,exclusion_contraints{ff}));
    end
    
    %
    empty_inclusion_constraints = cellfun(@(x)isempty(x),inclusion_constraints);
    for ff = 1 : n_features
        if empty_inclusion_constraints(ff)
            feature_flags(:,ff) = true;
        else
            feature_flags(:,ff) = ismember(feature_values(:,ff),inclusion_constraints{ff});
        end
    end
    intersection_flags = all(feature_flags,2);

    %
    values2combine = cell(n_features,1);
    for ff = 1 : n_features
        if empty_inclusion_constraints(ff)
            values2combine{ff} = nan;
        else
            values2combine{ff} = unique(feature_values(intersection_flags,ff));
        end
    end
    value_combinations = cell(n_features,1);
    [value_combinations{:}] = ndgrid(values2combine{:});
    value_combinations = cellfun(@(x)x(:),value_combinations,...
        'uniformoutput',false);

    %
    idcs2combine = cell(n_features,1);
    for ff = 1 : n_features
        if empty_inclusion_constraints(ff)
            idcs2combine{ff} = nan;
        else
            idcs2combine{ff} = find(ismember(feature_sets{ff},values2combine{ff}));
        end
    end
    idx_combinations = cell(n_features,1);
    [idx_combinations{:}] = ndgrid(idcs2combine{:});
    idx_combinations = cellfun(@(x)x(:),idx_combinations,...
        'uniformoutput',false);

    %
    conditions.n = max(cellfun(@(x)numel(x),value_combinations));
    for ff = 1 : n_features
        if empty_inclusion_constraints(ff)
            conditions.idcs.(feature_names{ff}) = ...
                repmat(1:numel(feature_sets{ff}),conditions.n,1);
            conditions.values.(feature_names{ff}) = ...
                repmat(feature_sets{ff},conditions.n,1);
        else
            conditions.idcs.(feature_names{ff}) = idx_combinations{ff};
            conditions.values.(feature_names{ff}) = value_combinations{ff};
        end
    end
    
    %
    conditions.idcs = struct2table(conditions.idcs);
    conditions.values = struct2table(conditions.values);
    
    %
    conditions2keep_flags = false(conditions.n,1);
    for kk = 1 : conditions.n
        feature_flags = false(n_trials,n_features);
        for ff = 1 : n_features
            feature_flags(:,ff) = ismember(...
                feature_values(:,ff),conditions.values(kk,:).(feature_names{ff}));
        end
        intersection_flags = all(feature_flags,2);
        conditions2keep_flags(kk) = sum(intersection_flags) > 0;
    end
    conditions.idcs = conditions.idcs(conditions2keep_flags,:);
    conditions.values = conditions.values(conditions2keep_flags,:);
    conditions.n = sum(conditions2keep_flags);
    
    %
    conditions.features.n = n_features;
    conditions.features.labels = feature_names;
end

