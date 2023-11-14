function [P_tX,P_Xt,pthat,features] = naivebayestimedecoder(X,opt)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here

    %% input parsing
    [n_timepoints,n_features,~] = size(X);

    %% output preallocation
    features = cell(n_features,1);
    for ff = 1 : n_features
        features{ff} = struct();
    end
    features = cell2mat(features);
    X_mus = nan(n_timepoints,n_features);
    P_Xt = nan(n_timepoints,n_features,opt.n_xpoints);
    P_tX = nan(n_timepoints,n_timepoints,opt.test.n_trials);
    pthat.mode = nan(n_timepoints,opt.test.n_trials);
    pthat.median = nan(n_timepoints,opt.test.n_trials);
    pthat.mean = nan(n_timepoints,opt.test.n_trials);

    %% overwrite MATLAB's builtin definition for the poisson PDF
    poisspdf = @(lambda,k) ...
        exp(k .* log(lambda + realmin) - lambda - gammaln(k + 1));

    %% estimate feature supports

    % iterate through features
    for ff = 1 : n_features
        if opt.verbose
            progressreport(ff,n_features,'estimating feature supports');
        end

        % parse feature
        x = squeeze(X(:,ff,:));

        % estimate feature span
        x_bounds = quantile(x(:),[0,1]+[1,-1]*.001).*(1+[-1,1]*.025);
        [~,x_edges] = histcounts(x(x>=x_bounds(1)&x<=x_bounds(2)),opt.n_xpoints);

        % update feature
        features(ff).idx = ff;
        features(ff).x_bounds = x_bounds;
        features(ff).x_edges = x_edges;
        features(ff).x_bw = range(x_bounds) / 10;
    end
    
    % concatenate feature supports
    X_edges = vertcat(features.x_edges);
    
    %% construct encoding models
    
    % iterate through features
    for ff = 1 : n_features
        if opt.verbose
            progressreport(ff,n_features,'constructing encoding models');
        end
        
        % parse feature
        x = squeeze(X(:,ff,:));
        x_bounds = features(ff).x_bounds;
        x_edges = features(ff).x_edges;
        x_bw = features(ff).x_bw;
        
        % compute tuning function
        X_mus(:,ff) = nanmean(x(:,opt.train.trial_idcs),2);
        
        % kernel definition
        x_kernel = normpdf(x_edges,mean(x_bounds),x_bw);
        x_kernel = x_kernel / nansum(x_kernel);
        
        % compute joint distribution
        if opt.assumepoissonmdl
            
            % store theoretical joint distribution
            P_Xt(:,ff,:) = poisspdf(X_mus(:,ff),x_edges(1:end-1));
        
            % temporal smoothing (acts as regularization)
%             P_Xt(:,ff,:) = nanconv2(squeeze(P_Xt(:,ff,:)),opt.t_kernel.pdf,1);
        else
            
            % preallocation
            X_counts = nan(n_timepoints,opt.train.n_trials,opt.n_xpoints);
            
            % iterate through training trials
            for kk = 1 : opt.train.n_trials
                train_idx = opt.train.trial_idcs(kk);
                
                % compute 2D histogram
                x_counts = histcounts2(1:n_timepoints,x(:,train_idx)',...
                    'xbinedges',1:n_timepoints+1,...
                    'ybinedges',x_edges);
                
                % re-nan what was nan before computing the 2D histogram
                nan_flags = isnan(x(:,train_idx));
                x_counts(nan_flags,:) = nan;
                
                % store current trial
                X_counts(:,kk,:) = x_counts;
            end
            
            % store average empirical joint distribution
            P_Xt(:,ff,:) = nanmean(X_counts,2);
            
            % 2D smoothing (acts as regularization)
            P_Xt(:,ff,:) = nanconv2(squeeze(P_Xt(:,ff,:)),x_kernel,x_kernel);
        end
        
%         if ff == 1
%             figure('position',[1.8000 41.8000 1.0224e+03 1.0288e+03]);
%         end
%         
%         subplot(3,1,1);
%         plot(opt.time,X_mus(:,ff));
%         
% %         subplot(3,1,2);
%         figure;
%         imagesc(squeeze(P_Xt(:,ff,:))');
%         set(gca,'ydir','normal');
        
        % zero fix (to prevent -inf issues when "logging" afterwards)
        if any(P_Xt(:,ff,:) == 0,'all')
            P_Xt(:,ff,:) = P_Xt(:,ff,:) + realmin;
        end

%         subplot(3,1,3);
%         imagesc(isnan(squeeze(P_Xt(:,ff,:)))');
%         pause(.1);
%         if any(isnan(P_Xt(:,ff,:)),'all')
%             a=1
%         end
        
        % update feature
        features(ff).x_mu = X_mus(:,ff);
        features(ff).p_Xc = squeeze(P_Xt(:,ff,:));
    end
    
    % for numerical reasons
    log_P_Xt = log(P_Xt);
    
    %% prior definition
    p_t = ones(n_timepoints,1) / n_timepoints;
    
    % for numerical reasons
    log_p_t = log(p_t);
    
    %% construct posteriors
    
    % iterate through test trials
    for kk = 1 : opt.test.n_trials
        if opt.verbose
            progressreport(kk,opt.test.n_trials,'constructing posteriors');
        end
        test_idx = opt.test.trial_idcs(kk);
        
        % iterate through time for the current test trial
        for tt = 1 : n_timepoints

            % fetch current observations
            x = X(tt,:,test_idx)';
            if all(isnan(x))
                continue;
            end
            
            % compute posterior for the current time point
            p_tX = decode(...
                x,X_edges,log_P_Xt,log_p_t,n_features,n_timepoints);

            % store posterior
            P_tX(tt,:,kk) = p_tX;
        end
        
        % flag valid time for the current trial
        test_time_flags = ~all(isnan(P_tX(:,:,kk)),2);

        % fetch single trial posteriors to compute point estimates
        P_tX_kk = P_tX(test_time_flags,:,kk);

        % posterior mode (aka MAP)
        [~,mode_idcs] = max(P_tX_kk,[],2);
        pthat.mode(test_time_flags,kk) = opt.time(mode_idcs);

        % posterior median
        median_flags = [false(sum(test_time_flags),1),...
            diff(cumsum(P_tX_kk,2) > .5,1,2) == 1];
        [~,median_idcs] = max(median_flags,[],2);
        pthat.median(test_time_flags,kk) = opt.time(median_idcs);

        % posterior mean (aka COM)
        P_tX_kk(isnan(P_tX_kk)) = 0;
        pthat.mean(test_time_flags,kk) = opt.time * P_tX_kk';
    end
    
    %% compute performance statistics
    
    % iterate through test trials
    for kk = 1 : opt.test.n_trials
        if opt.verbose
            progressreport(kk,opt.test.n_trials,'computing performance stats');
        end
        test_idx = opt.test.trial_idcs(kk);
    end
    
    % accuracy & precision stats..
end

function p_tX = decode(x,X_edges,log_P_Xt,log_p_t,n_features,n_timepoints)

    % index current observation
    [~,x_idcs] = min(abs(X_edges(:,1:end-1) - x),[],2);

    % preallocation
    log_p_tx = nan(n_features,n_timepoints);

    % iterate through features
    for ff = 1 : n_features

        % assume empirical encoding model
        log_p_tx(ff,:) = log_P_Xt(:,ff,x_idcs(ff));
    end
    
    % nan check
    nan_flags = all(isnan(log_p_tx),2) | isnan(x);
    if all(nan_flags)
        return;
    end

    % compute posterior by summing over log-likelihoods
    log_p_tX = log_p_t + sum(log_p_tx(~nan_flags,:))';
    
    % exponentiate to get back to probability
    p_tX = exp(log_p_tX - nanmax(log_p_tX));

    % normalization
    p_tX = p_tX / nansum(p_tX);
    
    if any(isnan(p_tX))
%         figure; plot(p_tX);
%         a=1
    end
end