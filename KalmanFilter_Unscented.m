function [ Mean, Var ] = KalmanFilter_Unscented( obs, init_state, init_var )
%KALMANFILTER Kalman filter a set of observations to give a set of means
% and covariances of the state

global Par;

% Get window length
L = size(obs, 1);

% If the length of the window is 0 then just return an empty array
if L == 0
    Mean = cell(0,1);
    Var = cell(0,1);
    return
end

% Create some arrays the usual KF quantities
PredMean = cell(L, 1);
PredVar = cell(L, 1);

Mean = cell(L, 1);
Var = cell(L, 1);

% Loop through time
for k = 1:L
    
    if k==1
        state = init_state;
        var = init_var;
    else
        state = Mean{k-1};
        var = Var{k-1};
    end
    
    % Unscented Transform
    [sig_pts_now, sig_wts, Np] = UnscentedTransform([state; zeros(4,1)], [var, zeros(4); zeros(4), Par.Q]);
    sig_pts_now(4,:) = max(sig_pts_now(4,:), Par.MinSpeed);
    
    % Project Forwards
    if Par.FLAG_DynMod == 0
        sig_pts_next = Par.A * sig_pts_now(1:4,:) + sig_pts_now(5:8,:);
    elseif Par.FLAG_DynMod == 1
        sig_pts_next = zeros(4,Np);
        for ii = 1:Np
            sig_pts_next(:,ii) = IntrinsicDynamicEvaluate(sig_pts_now(1:4,ii), sig_pts_now(5:8,ii));
        end
    end
    
    % Prediction mean and variance
    PredMean{k} = sig_pts_next * sig_wts;
	PredVar{k} = zeros(4);
    diff_next = bsxfun(@minus, sig_pts_next, PredMean{k});
    for ii = 1:Np
        PredVar{k} = PredVar{k} + sig_wts(ii) * diff_next(:,ii)*diff_next(:,ii)';
    end
    
    % Update step
    if isempty(obs{k})
        
        % No observation, so no update
        Mean{k} = PredMean{k};
        Var{k} = PredVar{k};
        
    else
        
        
        if Par.FLAG_ObsMod == 0
            
            % Normal KF update
            C = Par.C;
            y = obs{k} - C * PredMean{k};
            s = C * PredVar{k} * C' + Par.R;
            gain = PredVar{k} * C' / s;
            Mean{k} = PredMean{k} + gain * y;
            Var{k} = (eye(4)-gain*C) * PredVar{k};
            
        elseif Par.FLAG_ObsMod == 1
            
            % Send points though the observation model
            sig_pts_obs = zeros(2,Np);
            for ii = 1:Np
                [bng, rng] = cart2pol(sig_pts_next(1,ii), sig_pts_next(2,ii));
                sig_pts_obs(:,ii) = [bng; rng];
            end
            
            % Remaining stats
            obs_mean = sig_pts_obs * sig_wts;
            obs_var = Par.R;
            cross_covar = zeros(4,2);
            diff_obs = bsxfun(@minus, sig_pts_obs, obs_mean);
            for ii = 1:Np
                obs_var = obs_var + sig_wts(ii) * diff_obs(:,ii)*diff_obs(:,ii)';
                cross_covar = cross_covar + sig_wts(ii) * diff_next(:,ii)*diff_obs(:,ii)';
            end
            
            % Update
            gain = cross_covar / obs_var;
            Mean{k} = PredMean{k} + gain * (obs{k} - obs_mean);
            Var{k} = PredVar{k} - gain * obs_var * gain';
            
        end
        
    end
    
    Var{k} = 0.5*(Var{k}+Var{k}');
    
end

end