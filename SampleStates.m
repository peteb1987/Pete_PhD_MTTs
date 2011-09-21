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

if Par.FLAG_DynMod == 0
    [Mean, Var] = KalmanFilter(Obs, init_state, init_var);
elseif Par.FLAG_DynMod == 1
    [Mean, Var] = KalmanFilter_Unscented(Obs, init_state, init_var);
end

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
    
    if (k == L) && (t == Set.tracks(j).death-1)
        % Sampling last state in track and window
        mu = Mean{L};
        sigma = Var{L};
    else
        % Sampling state last in neither window or track 
        if (Par.FLAG_DynMod == 1)
            % Calculate A_{tt-1} and Q_{tt}
            
%             [A, Q] = IntrinsicDynamicLinearise(next_state);
%             sigma = inv(A' * (Q \ A) + inv(Var{k}));
%             mu = sigma * (A' * (Q \ next_state) + (Var{k} \ Mean{k}));

%             % generate sigma points
%             sig_pts = zeros(4,9);
%             sig_wts = zeros(9,1);
%             sig_wts(1) = -1/3;
%             for spi = 2:9
%                 col = mod(spi-2,4)+1;
%                 if spi<6, sgn=1; else sgn=-1; end
%                 sig_pts(:,spi) = sgn * Par.UQchol(:,col);
%                 sig_wts(spi) = 1/6;
%             end
%             
%             % Propagate points backwards
%             back_sig_pts = zeros(4,9);
%             for spi = 1:9
%                 back_sig_pts(:,spi) = IntrinsicDynamicInverse(next_state, sig_pts(:,spi));
%             end
%             back_mean = back_sig_pts * sig_wts;
%             back_less_mean = bsxfun(@minus, back_sig_pts, back_mean);
%             back_var = zeros(4,4);
%             for spi = 1:9
%                 back_var = back_var + sig_wts(spi) * back_less_mean(:,spi)*back_less_mean(:,spi)';
%             end
%             
%             sigma = inv(inv(back_var) + inv(Var{k}));
%             mu = sigma * (back_var\back_mean + Var{k}\Mean{k}); %#ok<MINV>

            % Unscented Transform
            [sig_pts_now, sig_wts, Np] = UnscentedTransform([Mean{k}; zeros(4,1)], [Var{k}, zeros(4); zeros(4), Par.Q]);
            sig_pts_now(4,:) = max(sig_pts_now(4,:), Par.MinSpeed);
            
            % Project Forwards
            for ii = 1:Np
                sig_pts_next(:,ii) = IntrinsicDynamicEvaluate(sig_pts_now(1:4,ii), sig_pts_now(5:8,ii));
            end
            
            % Calculate stats
            pred_mean = sig_pts_next * sig_wts;
            diff_next = bsxfun(@minus, sig_pts_next, pred_mean);
            diff_now = bsxfun(@minus, sig_pts_now(1:4,:), Mean{k});
            pred_covar = zeros(4,4);
            for ii = 1:Np
                pred_covar = pred_covar + sig_wts(ii) * diff_next(:,ii)*diff_next(:,ii)';
            end
            cross_covar = zeros(4);
            for ii = 1:Np
                cross_covar = cross_covar + sig_wts(ii) * diff_now(:,ii)*diff_next(:,ii)';
            end
            
            sigma = Var{k} - (cross_covar/pred_covar)*cross_covar';
            mu = Mean{k} + (cross_covar/pred_covar)*(next_state - pred_mean);
            
        elseif (Par.FLAG_DynMod == 0)
            sigma = inv(A' * (Q \ A) + inv(Var{k}));
            mu = sigma * (A' * (Q \ next_state) + (Var{k} \ Mean{k})); %#ok<MINV>
        end
    end
    
    sigma = (sigma+sigma')/2;
    
%     if t==2, figure(1); plot_gaussian_ellipsoid(Mean{L}(1:2), Var{L}(1:2,1:2)); end
    
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
