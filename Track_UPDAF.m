function [ Results ] = Track_UPDAF( detections, Observs, InitState )
%TRACK_PDAF Track targets using PDAF

global Par;
global Templates;

% Initialise arrays for results, intermediates and diagnostics
Results = cell(Par.T, 1);

if ~Par.FLAG_KnownInitStates
    % No knowledge of target starting positions
    InitEst = Templates.TrackSet;
    
else
    % Start with initial particle locations
    InitEst = Templates.TrackSet;
    InitEst.origin = ones(1, Par.NumTgts);
    InitEst.origin_time = ones(1, Par.NumTgts);
    for j = 1:Par.NumTgts
        if ~isempty(InitState{j})
            % Create a new track
            track = Templates.Track;
            track.birth = 0; track.death = 1; track.num = 1;
            track.state{1} = InitState{j};
            track.covar{1} = Par.KFInitVar*eye(4);
            track.assoc = 0;
            
            % Add it to the set
            InitEst.tracks = [InitEst.tracks; track];
            InitEst.N = InitEst.N + 1;
            InitEst.members = [InitEst.members; j];
        end
    end
end

% Loop through time
for t = 1:Par.T
    
    tic;
    
    disp('**************************************************************');
    disp(['*** Now processing frame ' num2str(t)]);
    
    if t==1
        [Results{t}] = UPDAFFrame(t, InitEst, Observs);
    else
        [Results{t}] = UPDAFFrame(t, Results{t-1}, Observs);
    end
    
    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');
    
end

end



function Est = UPDAFFrame(t, PrevEst, Observs)

global Par;

Est = PrevEst;

% Project tracks forward
for j = 1:Est.N
    if t == Est.tracks(j).death
        Est.tracks(j).death = Est.tracks(j).death + 1;
        Est.tracks(j).num = Est.tracks(j).num + 1;
        Est.tracks(j).state = [Est.tracks(j).state; {zeros(4,1)}];
        Est.tracks(j).smooth = [Est.tracks(j).smooth; {zeros(4,1)}];
        Est.tracks(j).covar = [Est.tracks(j).covar; {zeros(4,4)}];
        Est.tracks(j).assoc = [Est.tracks(j).assoc; 0];
    end
end

for j = 1:Par.NumTgts
    
    prev_state = Est.tracks(j).state{t-1 -Est.tracks(j).birth+1};
    prev_covar = Est.tracks(j).covar{t-1 -Est.tracks(j).birth+1};
    
    % Form sigma points of augmented state
    na = 8;
    kappa = 0;                     % IF YOU GET NON-POSITIVE DEFINITE COVARIANCES, INCREASE THIS (PD assured if k>0)
    sig_pts = zeros(na,2*na+1);
    sig_wts = zeros(2*na+1,1);
    sig_pts(1:4,1) = prev_state;
    sig_wts(1) = kappa/(na+kappa);
    mat_sq_rt = chol( (na+kappa)* [prev_covar, zeros(4); zeros(4), Par.Q] )';
    for spi = 2:2*na+1
        col = mod(spi-2,na)+1;
        if spi<na+2, sgn=1; else sgn=-1; end
        sig_pts(:,spi) = [prev_state; zeros(4,1)] + sgn * mat_sq_rt(:,col);
        sig_wts(spi) = 1/(2*(na+kappa));
    end
            
    % Propagate points forwards
    forward_sig_pts = zeros(4,17);
    if Par.FLAG_DynMod == 0
        forward_sig_pts = Par.A * sig_pts(1:4,:) + sig_pts(5:8,:);
    elseif Par.FLAG_DynMod == 1
        for spi = 1:17
            forward_sig_pts(:,spi) = IntrinsicDynamicEvaluate(sig_pts(1:4,spi), sig_pts(5:8,spi));
        end
    end
    
    % Calculate predictive state mean and variance
    pred_state = forward_sig_pts * sig_wts;
    state_diff = bsxfun(@minus, forward_sig_pts, pred_state);
    pred_covar = zeros(4,4);
    for spi = 1:17
        pred_covar = pred_covar + sig_wts(spi) * state_diff(:,spi)*state_diff(:,spi)';
    end
    
    % Propagate sigma points through observation function
    if Par.FLAG_ObsMod == 0
        obs_sig_pts = Par.C * forward_sig_pts;
    elseif Par.FLAG_ObsMod == 1
        [bng, rng] = cart2pol(forward_sig_pts(1, :), forward_sig_pts(2, :));
        obs_sig_pts = [bng; rng];
    end
    
    % Calculate predictive observation mean and variance
    mu = obs_sig_pts * sig_wts;
    innov = bsxfun(@minus, obs_sig_pts, mu);
    S = Par.R;
    for spi = 1:17
        S = S + sig_wts(spi) * innov(:,spi)*innov(:,spi)';
    end
    
    % Calculate cross correlation
    P = zeros(4,2);
    for spi = 1:17
        P = P + sig_wts(spi) * ( forward_sig_pts(:,spi)-pred_state )*innov(:,spi)';
    end
    
    % Gate
    ind = 1:Observs(t).N;
    innov = bsxfun(@minus, Observs(t).r(ind,:), mu')';
    if Par.FLAG_ObsMod == 1
        wrap_around1 = innov(1,:)>pi; innov(1, wrap_around1) = innov(1, wrap_around1) - 2*pi;
        wrap_around2 = innov(1,:)<-pi; innov(1, wrap_around2) = innov(1, wrap_around2) + 2*pi;
    end
    thresh1 = Par.GateSDs*sqrt(S(1,1)); thresh2 = Par.GateSDs*sqrt(S(2,2));
    test1 = abs(innov(1, :)) < thresh1; test2 = abs(innov(2, :)) < thresh2;
    indexes = find(test1&test2);
    test = false(1, length(ind));
    for i = indexes
        test(i) = ((innov(:,i)'/S)*innov(:,i) < Par.GateSDs^2);
    end
    validated = ind(test);
    N_gated = length(validated);
    gated_innov = innov(:, validated);
    
    % Calculate association posteriors and thus total innovation
	assocs = Par.PDetect * mvnpdf(gated_innov', [0 0], S);
    assoc_clut = (1-Par.PDetect)*Par.ExpClutObs*Par.ClutDens;
    norm_const = sum(assocs) + assoc_clut;
    assocs = assocs / norm_const;
    assoc_clut = assoc_clut / norm_const;
    tot_innov = gated_innov * assocs;
    
    % Calculate mean estimate
    W = P / S;
    est_mean = pred_state + W * tot_innov;
    
    % Calculate covariance estimate
    weighted_variance = zeros(2, 2, N_gated);
    for i = 1:N_gated
        weighted_variance(:,:,i) = assocs(i) * gated_innov(:, i) * gated_innov(:, i)';
    end
    middle_bit = sum(weighted_variance, 3) - tot_innov*tot_innov';
    est_covar = assoc_clut*pred_covar + (1-assoc_clut)*(pred_covar-W*S*W') + W*middle_bit*W';
    
    % Update track
    Est.tracks(j).state{t -Est.tracks(j).birth+1} = est_mean;
    Est.tracks(j).covar{t -Est.tracks(j).birth+1} = est_covar;

end

end

