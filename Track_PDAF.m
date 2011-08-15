function [ Results ] = Track_PDAF( detections, Observs, InitState )
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
        [Results{t}] = PDAFFrame(t, InitEst, Observs);
    else
        [Results{t}] = PDAFFrame(t, Results{t-1}, Observs);
    end
    
    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');
    
end

end



function Est = PDAFFrame(t, PrevEst, Observs)

global Par;

Est = PrevEst;

% Project tracks forward
for j = 1:Est.N
    if t == Est.tracks(j).death
        state = Par.A * Est.tracks(j).state{t-1-Est.tracks(j).birth+1};
        covar = Par.A * Est.tracks(j).covar{t-1-Est.tracks(j).birth+1} * Par.A' + Par.Q;
        Est.tracks(j).death = Est.tracks(j).death + 1;
        Est.tracks(j).num = Est.tracks(j).num + 1;
        Est.tracks(j).state = [Est.tracks(j).state; {state}];
        Est.tracks(j).covar = [Est.tracks(j).covar; {covar}];
        Est.tracks(j).assoc = [Est.tracks(j).assoc; 0];
    end
end

for j = 1:Par.NumTgts

    pred_state = Est.tracks(j).state{t -Est.tracks(j).birth+1};
    pred_covar = Est.tracks(j).covar{t -Est.tracks(j).birth+1};
    
    % Predict forwards mean and covariance
    if Par.FLAG_ObsMod == 0
        C = Par.C;
        mu = pred_state(1:2);
        
    elseif Par.FLAG_ObsMod == 1
        C = zeros(2, 4);
        x1 = pred_state(1); x2 = pred_state(2);
        C(1,1) = -x2/(x1^2+x2^2);
        C(1,2) = x1/(x1^2+x2^2);
        C(2,1) = x1/sqrt(x1^2+x2^2);
        C(2,2) = x2/sqrt(x1^2+x2^2);
        
        [bng, rng] = cart2pol(pred_state(1), pred_state(2));
        mu = [bng; rng];
        
    end
    
    S = C * pred_covar * C' + Par.R;
    S = 0.5*(S+S');

    % Gate
    ind = 1:Observs(t).N;
    innov = bsxfun(@minus, Observs(t).r(ind,:), mu')';
    if Par.FLAG_ObsMod == 1
        wrap_around1 = innov(1,:)>pi; innov(1, wrap_around1) = innov(1, wrap_around1) - 2*pi;
        wrap_around2 = innov(1,:)<-pi; innov(1, wrap_around2) = innov(1, wrap_around2) + 2*pi;
    end
    thresh1 = 5*sqrt(S(1,1)); thresh2 = 5*sqrt(S(2,2));
    test1 = abs(innov(1, :)) < thresh1; test2 = abs(innov(2, :)) < thresh2;
    indexes = find(test1&test2);
    test = false(1, length(ind));
    for i = indexes
        test(i) = ((innov(:,i)'/S)*innov(:,i) < 25);
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
    W = pred_covar * C' / S;
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

