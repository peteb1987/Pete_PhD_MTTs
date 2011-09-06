function [ Mean, Var ] = KalmanFilter( obs, init_state, init_var )
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

P = Par.P;

% Set C depending on model
if Par.FLAG_ObsMod == 0
    C = Par.C;
elseif Par.FLAG_ObsMod == 1
    C = zeros(2, 4);
end

% Set A depending on model
if Par.FLAG_DynMod == 0
    A = Par.A;
    Q = Par.Q;
elseif Par.FLAG_DynMod == 1
    A = zeros(4,4);
    Q = zeros(4,4);
end

% Loop through time
for k = 1:L
    
    if k==1
        state = init_state;
        var = init_var;
    else
        state = Mean{k-1};
        var = Var{k-1};
    end
    
    % Dynamic model linearisation
    if (Par.FLAG_DynMod == 1)
        [A, Q] = IntrinsicDynamicLinearise(state);
    end
    
    % Prediction step
    if Par.FLAG_DynMod == 0
        PredMean{k} = A * state;
    elseif Par.FLAG_DynMod == 1
        PredMean{k} = IntrinsicDynamicPredict(state);
    end
    PredVar{k} = A * var * A' + Q;
    
    % Update step
    if ~isempty(obs{k})
        
        % Observation associated with target
        if (Par.FLAG_ObsMod == 1)
            % Observation model linearisation
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

end