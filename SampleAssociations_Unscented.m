function [assoc, assoc_prob] = SampleAssociations_Unscented(j, t, L, Set, Observs, no_samp)
%SAMPLEASSOCIATIONS Sample current associations using unscented method

global Par;

% Initialise arrays
assoc = zeros(L, 1);
frame_prob = zeros(L, 1);
    
% Get t-L state and var
init_state = Set.tracks(j).state{t-L-Set.tracks(j).birth+1};
if ~Par.FLAG_RB
    init_var = Par.KFInitVar*eye(4);
else
    init_var = Set.tracks(j).covar{t-L-Set.tracks(j).birth+1};
end

% Set offset index
d = inf;                                                                   % d is the number of frames after the current one (tt) in which an observation is detected
for tt = t+1:Set.tracks(j).death-1
    ass = Set.tracks(j).assoc(tt-Set.tracks(j).birth+1);
    if ass ~= 0
        % Found a future observation - this means we're sampling a bridge
        d = tt - t;
        ass_next = Set.tracks(j).assoc(tt -Set.tracks(j).birth+1);
        y_next = Observs(tt).r(ass_next, :)';
    end
end

if ~isinf(d)
    Lpd = L + d;
else
    Lpd = L;
end

% Create a grand augmented noise covariance matrix
grand_mean = init_state;
grand_covar = init_var;
for k = 1:Lpd
    grand_covar = [grand_covar, zeros(4*k,4); zeros(4,4*k), Par.Q];
    grand_mean = [grand_mean; zeros(4,1)];
end

% Create Sigma Points
[ grand_sig_pts, sig_wts, Np ] = UnscentedTransform(grand_mean, grand_covar);

% Send points through motion model
XSigPts = cell(Lpd,1);
curr_pts = grand_sig_pts(1:4,:);
for k = 1:Lpd
    if Par.FLAG_DynMod == 0
        XSigPts{k} = Par.A * curr_pts + grand_sig_pts( 4*(k-1)+1:4*k ,:);
    elseif Par.FLAG_DynMod == 1
        XSigPts{k} = zeros(4,Np);
        for ii = 1:Np
            XSigPts{k}(:,ii) = IntrinsicDynamicEvaluate(curr_pts(:,ii), grand_sig_pts( 4*k+1:4*(k+1) ,ii));
        end
    end
    curr_pts = XSigPts{k};
end

% Send points though the observation model
YSigPts = cell(Lpd,1);
for k = 1:Lpd
    if Par.FLAG_ObsMod == 0
        YSigPts{k} = Par.C * XSigPts{k};
    elseif Par.FLAG_ObsMod == 1
        YSigPts{k} = zeros(2,Np);
        for ii = 1:Np
            [bng, rng] = cart2pol(XSigPts{k}(1,ii), XSigPts{k}(2,ii));
            YSigPts{k}(:,ii) = [bng; rng];
        end
    end
end

% Calculate the mean and covariance stats
XPredMeans = cell(Lpd,1);
YPredMeans = cell(Lpd,1);
for k = 1:Lpd
    XPredMeans{k} = XSigPts{k} * sig_wts;
    YPredMeans{k} = YSigPts{k} * sig_wts;
end

% Find last state of track (i.e. current time or when it dies)
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
    
    % Get required means
    p_y_now = YPredMeans{k};
    if ~isinf(d)
        p_y_next = YPredMeans{k+d};
    end
    
    % Calculate required covariances
    y_var_now = Par.R;
    diff_now = bsxfun(@minus, YSigPts{k}, YPredMeans{k});
    for ii = 1:Np
        y_var_now = y_var_now + sig_wts(ii) * diff_now(:,ii)*diff_now(:,ii)';
    end
    if ~isinf(d)
        y_var_next = Par.R;
        y_covar_future = zeros(2);
        diff_next = bsxfun(@minus, YSigPts{k+d}, YPredMeans{k+d});
        for ii = 1:Np
            y_var_next = y_var_next + sig_wts(ii) * diff_next(:,ii)*diff_next(:,ii)';
            y_covar_future = y_covar_future + sig_wts(ii) * diff_now(:,ii)*diff_next(:,ii)';
        end
    end
    
    % Calculate mean and variance
    if isinf(d)
        m = p_y_now;
        S = y_var_now;
    else
        m = p_y_now + (y_covar_future / y_var_next)*(y_next - p_y_next);
        S = y_var_now - (y_covar_future / y_var_next)*y_covar_future';
    end
    
    S = (S+S')/2;
    
    thresh1 = Par.GateSDs*sqrt(S(1,1));
    thresh2 = Par.GateSDs*sqrt(S(2,2));
    
    % Precalculate some things to speed up the next loop, which is over
    % 100's of observations within a loop over 100's of particles (i.e.
    % the bottleneck of the whole program)
    ind = 1:Observs(tt).N;
    innov = bsxfun(@minus, Observs(tt).r(ind,:), m')';
    if Par.FLAG_ObsMod == 1
        wrap_around1 = innov(1,:)>pi; innov(1, wrap_around1) = innov(1, wrap_around1) - 2*pi;
        wrap_around2 = innov(1,:)<-pi; innov(1, wrap_around2) = innov(1, wrap_around2) + 2*pi;
    end
    test1 = abs(innov(1, :)) < thresh1;
    test2 = abs(innov(2, :)) < thresh2;
    indexes = find(test1&test2);
    test = false(1, length(ind));
    for i = indexes
        test(i) = ((innov(:,i)'/S)*innov(:,i) < Par.GateSDs^2);
    end
    validated = ind(test);
    
    % Calculate weights
    ppsl_weights = zeros(N+1, 1);
    for i = validated
        ppsl_weights(i) = (Par.PDetect) * mvnpdf(Observs(tt).r(i, :), m', S);
    end
    
    %         %%%
    %         ppsl_weights = ppsl_weights * Par.PDetect;
    %
    %         ppsl_weights(N+1) = Par.ExpClutObs * Par.ClutDens * (1-Par.PDetect);
    %         %%%
    
    % Clutter
    ppsl_weights(N+1) = Par.ClutDens * (1-Par.PDetect);
    
%     if Par.FLAG_AlgType ~= 4
%         % Remove used ones
%         ppsl_weights(obs_assigned) = 0;
%     end
    
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
    
    % Did we find an observation?
    if ass==0
        d=d+1;
    else
        d=1;
        y_next = Observs(tt).r(ass, :)';
    end
    
end

% assoc_prob = sum(frame_prob) + init_state_prob;
assoc_prob = sum(frame_prob);

end

