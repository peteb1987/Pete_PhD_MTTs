function [assoc, assoc_prob] = SampleAssociations(j, t, L, Set, Observs, no_samp)
%SAMPLEASSOCIATIONS Sample current associations

global Par;

% Initialise arrays
assoc = zeros(L, 1);
frame_prob = zeros(L, 1);

% Set some local matrixes
R = Par.R;
Q = Par.Q;
invR = inv(Par.R);
invQ = inv(Par.Q);
A = Par.A;

% Get t-L state
x = Set.tracks(j).state{t-L-Set.tracks(j).birth+1};

% Set offset index
d = inf;                    % d is the number of frames after the current one (tt) in which an observation is detected
for tt = t+1:Set.tracks(j).death-1
    k = tt - (t-L);
    ass = Set.tracks(j).assoc(tt-Set.tracks(j).birth+1);
    if ass ~= 0
        % Found a future observation - this means we're sampling a bridge
        d = tt - t;
        
        % Calculate required quantities for this point
        next_y = Observs(tt).r(ass, :)';
        if Par.FLAG_ObsMod == 0
            
        elseif Par.FLAG_ObsMod == 1
            p_x = (A^k) * x;
            [p_bng, p_rng] = cart2pol(p_x(1), p_x(2));
            p_rngsq = p_rng^2;
            J = [-p_x(2)/p_rngsq, p_x(1)/p_rngsq, 0, 0; p_x(1)/p_rng, p_x(2)/p_rng, 0, 0];
            next_p_pol = [p_bng; p_rng];
            next_p_x = p_x;
            next_J = J;
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
    
    % Calculate deterministic prediction and jacobian
    p_x = (A^k) * x;
    if Par.FLAG_ObsMod == 0
        C = [1 0 0 0; 0 1 0 0];
    elseif Par.FLAG_ObsMod == 1
        [p_bng, p_rng] = cart2pol(p_x(1), p_x(2));
        p_rngsq = p_rng^2;
        J = [-p_x(2)/p_rngsq, p_x(1)/p_rngsq, 0, 0; p_x(1)/p_rng, p_x(2)/p_rng, 0, 0];
    end
    
    % Calculate the mean and variance
    if Par.FLAG_ObsMod == 0
        if d>L
            mu = C*p_x;
            S = R;
            for kk = 0:k-1
                S = S + C*(A^kk)*Q*(A^kk)'*C';
            end
        else
            xV = zeros(4);
            for kk = 0:k-1
                xV = xV + (A^kk)*Q*(A^kk)';
            end
            yV = R + C*xV*C';
            invSig = C'*(R\C) + ((A^d)'*C'/yV)*C*(A^d) + inv(xV);
            invS = invR - ((R\C)/invSig)*C'/R;
            S = inv(invS);
            mu = invS \ ( (((R\C)*(invSig\(A^d)')*C'/yV)*next_y) + (((R\(C/invSig))/xV)*(A^k)*x) );
        end
    elseif Par.FLAG_ObsMod == 1
        if d>L
            mu = [p_bng; p_rng];
            S = R;
            for kk = 0:k-1
                S = S + J*(A^k)*Q*(A^k)'*J';
            end
        else
            xV = zeros(4);
            for kk = 0:k-1
                xV = xV + (A^kk)*Q*(A^kk)';
            end
            yV = R + next_J * xV * next_J';
            invSig = J'*(R\J) + ((A^d)'*next_J'/yV)*next_J*(A^d) + inv(xV);
            invS = invR - ((R\J)/invSig)*J'/R;
            S = inv(invS);
            mu = [p_bng; p_rng] - J*p_x + invS \ ( (((R\J)*(invSig\(A^d)')*next_J'/yV)*(next_y-next_p_pol+next_J*next_p_x)) + (((R\(J/invSig))/xV)*(A^k)*x) );
        end
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
    
    % Remove used ones
    ppsl_weights(obs_assigned) = 0;
    
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
        next_y = Observs(tt).r(ass, :)';
        d=1;
        if Par.FLAG_ObsMod == 0
            
        elseif Par.FLAG_ObsMod == 1
            next_p_pol = [p_bng; p_rng];
            next_p_x = p_x;
            next_J = J;
        end
    end
    
end

assoc_prob = sum(frame_prob);

end

