function [assoc, assoc_prob] = SampleAssociations(j, t, L, Set, Observs, no_samp)
%SAMPLEASSOCIATIONS Sample current associations

global Par;

% Initialise arrays
assoc = zeros(L, 1);
frame_prob = zeros(L, 1);

% Set some local matrixes
R = Par.R;
invR = inv(Par.R);
if Par.FLAG_DynMod == 0
    A = Par.A;
    Q = Par.Q;
elseif Par.FLAG_DynMod == 1
    A = zeros(4);
    Q = zeros(4);
end
    
% Get t-L state and var
if ~Par.FLAG_RB
    x = Set.tracks(j).state{t-L-Set.tracks(j).birth+1};
%     init_state_prob = 0;
    init_var = zeros(4);
else
    x = Set.tracks(j).state{t-L-Set.tracks(j).birth+1};
    init_var = Set.tracks(j).covar{t-L-Set.tracks(j).birth+1};
%     x = mvnrnd(Set.tracks(j).state{t-L-Set.tracks(j).birth+1}', Set.tracks(j).covar{t-L-Set.tracks(j).birth+1})';
%     init_state_prob = log(mvnpdf(x', Set.tracks(j).state{t-L-Set.tracks(j).birth+1}', Set.tracks(j).covar{t-L-Set.tracks(j).birth+1})');
end

% Set offset index
d = inf;                    % d is the number of frames after the current one (tt) in which an observation is detected
for tt = t+1:Set.tracks(j).death-1
    k = tt - (t-L);
    ass = Set.tracks(j).assoc(tt-Set.tracks(j).birth+1);
    if ass ~= 0
        % Found a future observation - this means we're sampling a bridge
        d = tt - t;
        
        % Calculate required quantities for this point
        y_next = Observs(tt).r(ass, :)';
        if Par.FLAG_DynMod == 0
            p_x = (A^k) * x;
        else
            p_x = x;
            for kk=1:k
                p_x = IntrinsicDynamicPredict(p_x);
            end
        end
        p_x_next = p_x;
        
        if Par.FLAG_ObsMod == 0
            C_next = Par.C;
            p_y_next = C_next * p_x_next;
        elseif Par.FLAG_ObsMod == 1
            [p_bng, p_rng] = cart2pol(p_x(1), p_x(2));
            p_rngsq = p_rng^2;
            C_next = [-p_x(2)/p_rngsq, p_x(1)/p_rngsq, 0, 0; p_x(1)/p_rng, p_x(2)/p_rng, 0, 0];
            p_y_next = [p_bng; p_rng];
        end
        
    end
end

% Find last state of track (i.e. current time of when it dies)
last = min(t, Set.tracks(j).death - 1);

% Loop through time
for tt = last:-1:t-L+1
    
    k = tt - (t-L);
    N = Observs(tt).N;
    
    % Create a list of used associations
    obs_assigned = [];
    for jj = 1:Set.N
        if (j ~= jj) && ((Set.tracks(j).birth <= tt)&&(Set.tracks(j).death > tt))
            ass = Set.tracks(jj).assoc(tt-Set.tracks(jj).birth+1);
            if ass > 0
                obs_assigned = [obs_assigned, ass];
            end
        end
    end
    
    % Calculate deterministic prediction
    if Par.FLAG_DynMod == 0
        p_x = (A^k) * x;
    else
        p_x = x;
        for kk=1:k
            p_x = IntrinsicDynamicPredict(p_x);
        end
    end
    
    % Calculate observation model jacobian
    if Par.FLAG_ObsMod == 0
        C = [1 0 0 0; 0 1 0 0];
        p_y = C*p_x;
    elseif Par.FLAG_ObsMod == 1
        [p_bng, p_rng] = cart2pol(p_x(1), p_x(2));
        p_y = [p_bng; p_rng];
        p_rngsq = p_rng^2;
        C = [-p_x(2)/p_rngsq, p_x(1)/p_rngsq, 0, 0; p_x(1)/p_rng, p_x(2)/p_rng, 0, 0];
    end
    
    % At this stage we should have p_x, p_y, and C calculated, the
    % predicted state, observation and linear observation matrix. Also the
    % same three suffixed by _next which give the same values d steps
    % forward, and y_next.
    
    if Par.FLAG_DynMod == 0
        Q_before = (A^k)*init_var*(A^k)';
        for kk = 0:k-1
            Q_before = Q_before + (A^kk)*Q*(A^kk)';
        end
    elseif Par.FLAG_DynMod == 1
        [~, Q_before] = IntrinsicDynamicCompoundStats(k, x, init_var);
%         Q_before = Q_before + 1E-3 * eye(4);
    end
    
    % Calculate the mean and variance
    if d > L
        % No known later observations
        mu = p_y;
        S = R + C*Q_before*C';
        
    else
        % Known later observation
        if Par.FLAG_DynMod == 0
            A_after = A^d;
            Q_after = zeros(4);
            for dd = 0:d
                Q_after = Q_after + (A^dd)*Q*(A^dd)';
            end
        elseif Par.FLAG_DynMod == 1
            [A_after, Q_after] = IntrinsicDynamicCompoundStats(d, p_x, zeros(4));
        end
        R_after = R + C_next*Q_after*C_next';
        invSigma = C'*(R\C) + (A_after'*C_next'/R_after)*C_next*A_after + inv(Q_before);
        invS = invR - ((R\C)/invSigma)*C'/R;
        S = inv(invS);
        mu = p_y - C*p_x + (invS\(R\C)/invSigma)*( A_after'*(C_next'/R_after)*(y_next-p_y_next+C_next*p_x_next) + Q_before\p_x );
        
    end
    
    S = (S+S')/2;
    
    thresh1 = 5*sqrt(S(1,1));
    thresh2 = 5*sqrt(S(2,2));
    
    % Precalculate some things to speed up the next loop, which is over
    % 100's of observations within a loop over 100's of particles (i.e.
    % the bottleneck of the whole program)
    ind = 1:Observs(tt).N;
    innov = bsxfun(@minus, Observs(tt).r(ind,:), mu')';
    if Par.FLAG_ObsMod == 1
        wrap_around1 = innov(1,:)>pi; innov(1, wrap_around1) = innov(1, wrap_around1) - 2*pi;
        wrap_around2 = innov(1,:)<-pi; innov(1, wrap_around2) = innov(1, wrap_around2) + 2*pi;
    end
    test1 = abs(innov(1, :)) < thresh1;
    test2 = abs(innov(2, :)) < thresh2;
    indexes = find(test1&test2);
    test = false(1, length(ind));
    for i = indexes
        test(i) = ((innov(:,i)'/S)*innov(:,i) < 25);
    end
    validated = ind(test);
    
    % Calculate weights
    ppsl_weights = zeros(N+1, 1);
    for i = validated
        ppsl_weights(i) = (Par.PDetect) * mvnpdf(Observs(tt).r(i, :), mu', S);
    end
    
    %         %%%
    %         ppsl_weights = ppsl_weights * Par.PDetect;
    %
    %         ppsl_weights(N+1) = Par.ExpClutObs * Par.ClutDens * (1-Par.PDetect);
    %         %%%
    
    % Clutter
    ppsl_weights(N+1) = Par.ClutDens * (1-Par.PDetect);
    
    if Par.FLAG_AlgType ~= 4
        % Remove used ones
        ppsl_weights(obs_assigned) = 0;
    end
    
    % Normalise
    ppsl_weights = ppsl_weights/sum(ppsl_weights);
    
    % Set minimum value for clutter as 1-PDetect
    %             if (d>L)&&(ppsl_weights(N+1)<(1-Par.PDetect))
    if (ppsl_weights(N+1)<(1-Par.PDetect))
        ppsl_weights(1:N) = Par.PDetect*ppsl_weights(1:N)/sum(ppsl_weights(1:N));
        ppsl_weights(N+1) = 1-Par.PDetect;
    end
    
    if ~no_samp
        % Sample
        ass = randsample(N+1, 1, true, ppsl_weights);
    else
        % Get the current value from the track
        ass = Set.tracks(j).assoc(tt-Set.tracks(j).birth+1);
        if ass == 0
            ass = N+1;
        end
    end
    
    % Find probability
    frame_prob(k) = log(ppsl_weights(ass));
    
    if ass == N+1
        ass = 0;
    end
    
    assoc(k) = ass;
    
    % Store numbers for next step
    if ass==0
        d=d+1;
    else
        d=1;
        y_next = Observs(tt).r(ass, :)';
        p_y_next = p_y;
        p_x_next = p_x;
        C_next = C;
    end
    
end

% assoc_prob = sum(frame_prob) + init_state_prob;
assoc_prob = sum(frame_prob);

end

