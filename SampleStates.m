function [State, state_prob, Mean, Var] = SampleStates(j, t, L, Set, Observs, no_samp)
%SAMPLESTATES Sample current states

global Par;

% Only need examine those which are present after t-L
if Set.tracks(j).death < t-L+1
    state_prob = 0;
    State = [];
    Var = [];
    return
end

% Set A and Q depending on model
if Par.FLAG_DynMod == 0
    A = Par.A;
    Q = Par.Q;
elseif Par.FLAG_DynMod == 1
    A = zeros(4,4);
    Q = zeros(4,4);
end

% How long should the KF run for?
last = min(t, Set.tracks(j).death - 1);
first = max(t-L+1, Set.tracks(j).birth+1);
num = last - first + 1;

% Initialise an array of states
State = cell(num, 1);
frame_prob = zeros(L, 1);

% Draw up a list of associated hypotheses
Obs = ListAssocObservs(last, num, Set.tracks(j), Observs);

% Run a Kalman filter for the target
init_state = Set.tracks(j).state{first-1 -Set.tracks(j).birth+1};
if ~Par.FLAG_RB
    init_var = Par.KFInitVar*eye(4);
else
    init_var = Set.tracks(j).covar{first-1 -Set.tracks(j).birth+1};
end
[Mean, Var] = KalmanFilter(Obs, init_state, init_var);

if Par.FLAG_RB
    State = Mean;
    state_prob = 0;
    return
end

% Loop through time
for k = L:-1:1
    
    tt = t-L+k;
    
    if (k == L) && (t < Set.tracks(j).death-1)
        % Sampling last state in window, but not track (e.g. bridging move)
        next_state = Set.tracks(j).state{last+1 -Set.tracks(j).birth+1};
    elseif k<L
        % Sampling state last in neither window or track
        next_state = State{k+1};
    end
    
    if (Par.FLAG_DynMod == 1)
        % Calculate A_{tt-1} and Q_{tt}
        [A, Q] = IntrinsicDynamicLinearise(Mean{k});
    end
    
    if (k == L) && (t == Set.tracks(j).death-1)
        % Sampling last state in track and window
        mu = Mean{L};
        sigma = Var{L};
    else
        % Sampling state last in neither window or track
        sigma = inv(A' * (Q \ A) + inv(Var{k}));
        mu = sigma * (A' * (Q \ next_state) + (Var{k} \ Mean{k})); %#ok<MINV>
    end
    
    sigma = (sigma+sigma')/2;
    
    if ~no_samp
        % Sample
        State{k} = mvnrnd(mu', sigma)';
        
        if (Par.FLAG_DynMod == 1)
            % Limit speed
            State{k}(4) = max(State{k}(4), Par.MinSpeed);
        end
        
    else
        % Get current value from the track
        State{k} = Set.tracks(j).state{tt -Set.tracks(j).birth+1};
    end
    
    frame_prob(k) = log( mvnpdf(State{k}', mu', sigma) );
    
end

state_prob = sum(frame_prob);

end
