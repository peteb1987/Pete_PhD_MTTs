function [ BackMean, BackVar ] = KalmanSmoother( obs, init_state, init_var )
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

% Set C depending on model
if Par.FLAG_ObsMod == 0
    % Fixed
    C = Par.C;

elseif Par.FLAG_ObsMod == 1
    % Use EKF approximation, calculated at each step
    C = zeros(2, 4);
    
end

% Loop through time
for k = 1:L
    
    % Prediction step
    
    if k==1
        PredMean{1} = Par.A * init_state;
        PredVar{1} = Par.A * init_var * Par.A' + Par.Q;
    else
        PredMean{k} = Par.A * Mean{k-1};
        PredVar{k} = Par.A * Var{k-1} * Par.A' + Par.Q;
    end
    
    % Update step
    
    if ~isempty(obs{k})
        
        % Observation associated with target
        
        if (Par.FLAG_ObsMod == 1)
            
            % Linearisation
            x1 = PredMean{k}(1);
            x2 = PredMean{k}(2);
            C(1,1) = -x2/(x1^2+x2^2);
            C(1,2) = x1/(x1^2+x2^2);
            C(2,1) = x1/sqrt(x1^2+x2^2);
            C(2,2) = x2/sqrt(x1^2+x2^2);
            
        end
        
        % Innovation
        if Par.FLAG_ObsMod == 0
            y = obs{k} - C * PredMean{k};
        elseif Par.FLAG_ObsMod == 1
            [bng, rng] = cart2pol(PredMean{k}(1), PredMean{k}(2));
            y = obs{k} - [bng; rng];
            if y(1) > pi
                y(1) = y(1) - 2*pi;
            elseif y(1) < -pi
                y(1) = y(1) + 2*pi;
            end
        end
        
        % Remaining KF calculations
        s = C * PredVar{k} * C' + Par.R;
        gain = PredVar{k} * C' / s;
        Mean{k} = PredMean{k} + gain * y;
        Var{k} = (eye(4)-gain*C) * PredVar{k};
        
    else
        % No observation, so no update
        Mean{k} = PredMean{k};
        Var{k} = PredVar{k};
        
    end
    
    Var{k} = 0.5*(Var{k}+Var{k}');
    
end

BackMean = cell(L, 1);
BackVar = cell(L, 1);
BackMean{L} = Mean{L}; BackVar{L} = Var{L};
% Smoothing
for k = L-1:-1:1
    
    G = Var{k}*Par.A'/PredVar{k+1};
    BackMean{k} = Mean{k} + G * (BackMean{k+1}-PredMean{k+1});
    BackVar{k} = Var{k} + G * (BackVar{k+1}-PredVar{k+1}) * G';
    
end

end