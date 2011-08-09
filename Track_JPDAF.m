function [ Results ] = Track_JPDAF( detections, Observs, InitState )
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
        [Results{t}] = PDAFFrame(t, t, InitEst, Observs);
    else
        [Results{t}] = PDAFFrame(t, min(t,Par.L), Results{t-1}, Observs);
    end
    
    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');
    
end

end



function Est = PDAFFrame(t, L, PrevEst, Observs)

global Par;

Est = PrevEst;

% Project tracks forward
for j = 1:Est.N
    if t == Est.tracks(j).death
        state = Par.A * Est.tracks(j).state{t-1-Est.tracks(j).birth+1};
        covar = Par.A * Est.tracks(j).covar{t-1-Est.tracks(j).birth+1} * Par.A + Par.Q;
        Est.tracks(j).death = Est.tracks(j).death + 1;
        Est.tracks(j).num = Est.tracks(j).num + 1;
        Est.tracks(j).state = [Est.tracks(j).state; {state}];
        Est.tracks(j).covar = [Est.tracks(j).covar; {covar}];
        Est.tracks(j).assoc = [Est.tracks(j).assoc; 0];
    end
end

pred_state = cell(Par.NumTgts,1);
pred_covar = cell(Par.NumTgts,1);
C = cell(Par.NumTgts,1);
S = cell(Par.NumTgts,1);
gated_innov = cell(Par.NumTgts,1);
validated = cell(Par.NumTgts,1);

for j = 1:Par.NumTgts
    
    pred_state{j} = Est.tracks(j).state{t -Est.tracks(j).birth+1};
    pred_covar{j} = Est.tracks(j).covar{t -Est.tracks(j).birth+1};
    
    % Predict forwards mean and covariance
    if Par.FLAG_ObsMod == 0
        C{j} = Par.C;
        mu = pred_state{j}(1:2);
        
    elseif Par.FLAG_ObsMod == 1
        C{j} = zeros(2, 4);
        x1 = pred_state(1); x2 = pred_state(2);
        C{j}(1,1) = -x2/(x1^2+x2^2);
        C{j}(1,2) = x1/(x1^2+x2^2);
        C{j}(2,1) = x1/sqrt(x1^2+x2^2);
        C{j}(2,2) = x2/sqrt(x1^2+x2^2);
        
        [bng, rng] = cart2pol(pred_state(1), pred_state(2));
        mu = [bng; rng];
        
    end
    
    S{j} = C{j} * pred_covar{j} * C{j}' + Par.R;
    S{j} = 0.5*(S{j}+S{j}');
    
    % Gate
    ind = 1:Observs(t).N;
    innov = bsxfun(@minus, Observs(t).r(ind,:), mu')';
    if Par.FLAG_ObsMod == 1
        wrap_around1 = innov(1,:)>pi; innov(1, wrap_around1) = innov(1, wrap_around1) - 2*pi;
        wrap_around2 = innov(1,:)<-pi; innov(1, wrap_around2) = innov(1, wrap_around2) + 2*pi;
    end
    thresh1 = 5*sqrt(S{j}(1,1)); thresh2 = 5*sqrt(S{j}(2,2));
    test1 = abs(innov(1, :)) < thresh1; test2 = abs(innov(2, :)) < thresh2;
    indexes = find(test1&test2);
    test = false(1, length(ind));
    for i = indexes
        test(i) = ((innov(:,i)'/S{j})*innov(:,i) < 25);
    end
    validated{j} = ind(test);
    gated_innov{j} = innov(:, validated{j});
    
    % Calculate target-association likelihoods
    likes{j} = mvnpdf(gated_innov{j}', [0 0], S{j});
    
end

% Check whether targets are all independent
indep = true;
for j = 1:Par.NumTgts
    for jj = j+1:Par.NumTgts
        if ~isempty(intersect(validated{j}, validated{jj}))
            indep = false;
            break
        end
    end
    if indep == false
        break
    end
end

% Find joint association probabilities
if ~indep
    num_valid = cellfun(@length, validated)+1;
    assoc_list = zeros(prod(num_valid), Par.NumTgts);
    
    for j = 1:Par.NumTgts
        
        inner_rep = prod(num_valid(1:j-1));
        outer_rep = prod(num_valid(j+1:end));
        
        temp = repmat([0; validated{j}']',inner_rep,1);
        temp = temp(:);
        
        assoc_list(:,j) = repmat(temp, outer_rep, 1);
        
    end
    
    joint_assoc = zeros(prod(num_valid), 1);
    
    % loop through and calculate joint probs
    for ii = 1:prod(num_valid)
        
        % Check for association collisions
        repeats = repval(assoc_list(ii, :));
        repeats(repeats==0)=[];
        
        if isempty(repeats)
            joint_assoc(ii) = 1;
            for j=1:Par.NumTgts
                if assoc_list(ii,j)==0
                    joint_assoc(ii) = joint_assoc(ii) * Par.ClutDens;
                else
                    joint_assoc(ii) = joint_assoc(ii) * likes{j}(validated{j}==assoc_list(ii,j));
                end
            end
        else
%             joint_assoc(ii) = 0;
        end
        
    end
    
end

for j = 1:Par.NumTgts
    
    % Caluclate association probs
    if indep
        assocs = Par.PDetect * mvnpdf(gated_innov{j}', [0 0], S{j});
        assoc_clut = Par.ClutDens * (1-Par.PDetect) / Par.PDetect;
    else
        assocs = zeros(num_valid(j)-1,1);
        for ii = 1:num_valid(j)-1
            ass = validated{j}(ii);
            assocs(ii) = sum(joint_assoc(assoc_list(:,j)==ass));
        end
        assoc_clut = sum(joint_assoc(assoc_list(:,j)==0));
    end
    
    % Calculate total innovation
    norm_const = sum(assocs) + assoc_clut;
    assocs = assocs / norm_const;
    assoc_clut = assoc_clut / norm_const;
    tot_innov = gated_innov{j} * assocs;
    
    % Calculate mean estimate
    W = pred_covar{j} * C{j}' / S{j};
    est_mean = pred_state{j} + W * tot_innov;
    
    % Calculate covariance estimate
    weighted_variance = zeros(2, 2, length(validated{j}));
    for i = 1:length(validated{j})
        weighted_variance(:,:,i) = assocs(i) * gated_innov{j}(:, i)' * gated_innov{j}(:, i);
    end
    middle_bit = sum(weighted_variance, 3) - tot_innov'*tot_innov;
    est_covar = assoc_clut*pred_covar{j} + (1-assoc_clut)*(pred_covar{j}-W*S{j}*W') + W*middle_bit*W';
    
    % Update track
    Est.tracks(j).state{t -Est.tracks(j).birth+1} = est_mean;
    Est.tracks(j).covar{t -Est.tracks(j).birth+1} = est_covar;
    
end

end

