function [P_tX,P_Xt,pthat,features] = naivebayesdecoder_simpler(tensor,opt)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    %% input parsing
    [n_timepoints,n_features,n_trials] = size(tensor);
    stimulus_set = unique(opt.stimulus);
    condition_set = unique(opt.condition);
    n_conditions = numel(condition_set);
    
    %% output preallocation
    features = cell(n_features,1);
    for ff = 1 : n_features
        features{ff} = struct();
    end
    features = cell2mat(features);
    pthat.mode = nan(n_timepoints,opt.test.n_trials);
    pthat.median = nan(n_timepoints,opt.test.n_trials);
    pthat.mean = nan(n_timepoints,opt.test.n_trials);
    P_Xt = nan(n_timepoints,n_features,opt.n_xpoints);
    P_tX = nan(n_timepoints,n_timepoints,opt.test.n_trials);
    
    %% estimate feature spans

    % iterate through features
    for ff = 1 : n_features
        progressreport(ff,n_features,'estimating feature spans');

        % parse feature
        X = squeeze(tensor(:,ff,:));
        
        % estimate feature span
        x_bounds = quantile(X(:),[0,1]+[1,-1]*.01).*(1+[-1,1]*.05);
        [~,x_edges] = histcounts(X(X>=x_bounds(1)&X<=x_bounds(2)),opt.n_xpoints);
        
        % update feature
        features(ff).idx = ff;
        features(ff).x_bounds = x_bounds;
        features(ff).x_edges = x_edges;
        features(ff).x_bw = range(x_bounds) / 10;
    end
    
    %% compute condition-specific descriptive statistics

    % preallocation
    x_mus = nan(n_features,n_conditions);
    x_sigs = nan(n_features,n_conditions);
    
    % iterate through features
    for ff = 1 : n_features
        progressreport(ff,n_features,'computing condition-specific stats');
        
        % parse feature
        X = squeeze(tensor(:,ff,:));
        
        % iterate through conditions
        for cc = 1 : n_conditions
            cond_flags = opt.condition == condition_set(cc);

            % compute condition mean
            x_mus(ff,cc) = nanmean(X(:,cond_flags),[1,2]);
            
            % compute condition standard deviation
            x_sigs(ff,cc) = nanstd(X(:,cond_flags),0,[1,2]);
        end
        
        % center condition means on the training one
        x_mus(ff,:) = x_mus(ff,:) - ...
            nanmean(X(:,opt.train.trial_idcs),[1,2]);
        
        % normalize condition standard deviations by the training one
        x_sigs(ff,:) = x_sigs(ff,:) / ...
            nanstd(X(:,opt.train.trial_idcs),0,[1,2]);
    end
        
    %% construct encoding models

    % preallocation
    X_mus = nan(n_timepoints,n_features);
    
    % iterate through features
    for ff = 1 : n_features
        progressreport(ff,n_features,'constructing encoding models');
        
        % parse feature
        X = squeeze(tensor(:,ff,:));
        x_bounds = features(ff).x_bounds;
        x_edges = features(ff).x_edges;
        x_bw = features(ff).x_bw;
        
        % compute tuning function
        X_mus(:,ff) = nanmean(X(:,opt.train.trial_idcs),2);
        
        % kernel definition
        t_kernel = normpdf(linspace(0,1,n_timepoints),.5,.015);
        t_kernel = t_kernel / nansum(t_kernel);
        x_kernel = normpdf(x_edges,mean(x_bounds),x_bw);
        x_kernel = x_kernel / nansum(x_kernel);
        
        % preallocation
        p_Xt = nan(n_timepoints,opt.train.n_trials,opt.n_xpoints);

        % iterate through training trials
        for kk = 1 : opt.train.n_trials
            train_idx = opt.train.trial_idcs(kk);
            
            % compute likelihood
            x_counts = histcounts2(1:n_timepoints,X(:,train_idx)',...
                'xbinedges',1:n_timepoints+1,...
                'ybinedges',x_edges);
            p_Xt(:,kk,:) = conv2(1,x_kernel,x_counts,'same');
        end
        
        % store average joint distribution
        P_Xt(:,ff,:) = nanmean(p_Xt,2);
    end

    % normalization
    P_Xt = P_Xt ./ nansum(P_Xt,3);
    
    %% select features based on condition-specific statistics
    feature_flags = all(x_sigs > 0 & ~isnan(x_sigs) & ~isinf(x_sigs),2);
    n_features = sum(feature_flags);
    features = features(feature_flags);
    x_mus = x_mus(feature_flags,:);
    x_sigs = x_sigs(feature_flags,:);
    X_mus = X_mus(:,feature_flags);
    P_Xt = P_Xt(:,feature_flags,:);
    
    %% construct posteriors

    % overwrite poisson PDF definition
    poisspdf = @(lambda,k) exp(k .* log(lambda) - lambda - gammaln(k + 1));
    
    % prior definition
    p_t = ones(n_timepoints,1) / n_timepoints;
    
    % iterate through test trials
    for kk = 1 : opt.test.n_trials
        progressreport(kk,opt.test.n_trials,'constructing posteriors');
        test_idx = opt.test.trial_idcs(kk);
        test_time_idcs = 1 : floor(n_timepoints * ...
            opt.stimulus(test_idx) / stimulus_set(end));

        % iterate through true time for the current test trial
        for tt = test_time_idcs
            
            % fetch current observations
            x = tensor(tt,feature_flags,test_idx)';
            
            % compute likelihoods of the current observations
            if opt.assumepoissonmdl
                
                % assume a features are poisson-distributed
                p_tx = poisspdf(X_mus',round(x));
            else
                
                % index current observation
                x_edges = vertcat(features.x_edges);
                [~,x_idcs] = min(abs(x_edges(:,1:end-1) - x),[],2);
                
                % preallocation
                p_tx = nan(n_features,n_timepoints);
                
                % iterate through features
                for ff = 1 : n_features
                    
                    % assume empirical encoding model
                    p_tx(ff,:) = P_Xt(:,ff,x_idcs(ff));
                end
            end
            
            % normalization
            p_tx = p_tx ./ nansum(p_tx,2);
            nan_flags = all(isnan(p_tx),2);
            if all(nan_flags)
                continue;
            end
            
            % compute posterior (with numerical precision issues in mind)
            fudge = 1 + 1 / n_features;
            p_tX = p_t .* prod(p_tx(~nan_flags,:) * n_timepoints + fudge,1)';
            p_X = nansum(p_tX);
            P_tX(tt,:,kk) = p_tX / p_X;
        end
        
        % fetch single trial posteriors to compute point estimates
        P_tX_kk = P_tX(test_time_idcs,:,kk);
        
        % posterior mode (aka MAP)
        [~,mode_idcs] = max(P_tX_kk,[],2);
        pthat.mode(test_time_idcs,kk) = opt.time(mode_idcs);
        
        % posterior median
        median_flags = [false(test_time_idcs(end),1),...
            diff(cumsum(P_tX_kk,2) > .5,1,2) == 1];
        [~,median_idcs] = max(median_flags,[],2);
        pthat.median(test_time_idcs,kk) = opt.time(median_idcs);
        
        % posterior mean (aka COM)
        P_tX_kk(isnan(P_tX_kk)) = 0;
        pthat.mean(test_time_idcs,kk) = opt.time * P_tX_kk';
    end
end