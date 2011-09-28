function [ Results ] = Track_MCMC( detections, Observs, InitState )
%TRACK_MCMC Track targets using MCMC-PF

global Par;
global Templates;

% Initialise arrays for results, intermediates and diagnostics
Results = cell(Par.T, 1);
Chains = cell(Par.T, 1);
BestEsts = cell(Par.T, 1);
MoveTypes = cell(Par.T, 1);

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
            track.smooth{1} = InitState{j};
            track.covar{1} = Par.KFInitVar*eye(4);
            track.assoc = 0;
            
            % Add it to the set
            InitEst.tracks = [InitEst.tracks; track];
            InitEst.N = InitEst.N + 1;
            InitEst.members = [InitEst.members; j];
        end
    end
end

InitChain = struct( 'particles', cell(1), 'posteriors', zeros(1) );
InitChain.particles{1} = InitEst;

% Loop through time
for t = 1:Par.T
    
    tic;
    
    disp('**************************************************************');
    disp(['*** Now processing frame ' num2str(t)]);
    
    if t==1
        [Chains{t}, BestEsts{t}, MoveTypes{t}] = MCMCFrame(t, t, {InitChain}, InitEst, Observs);
    else
        [Chains{t}, BestEsts{t}, MoveTypes{t}] = MCMCFrame(t, min(t,Par.L), Chains(1:t), BestEsts{t-1}, Observs);
    end
    
    Results{t} = Chains{t};
    
    disp(['*** Correct associations at frame ' num2str(t-min(t,Par.L)+1) ': ' num2str(detections(t-min(t,Par.L)+1,:))]);
    assoc = [];
    for j = 1:Par.NumTgts
        get_ass = cellfun(@(x) x.tracks(j).assoc(t-min(t,Par.L)+1 -x.tracks(j).birth+1), Chains{t}.particles);
        mode_ass = mode(get_ass);
        assoc = [assoc, mode_ass];
    end
    disp(['*** Modal associations at frame ' num2str(t-min(t,Par.L)+1) ': ' num2str(assoc)]);
    assoc = arrayfun(@(x) x.assoc(t-min(t,Par.L)+1 -x.birth+1), BestEsts{t}.tracks)';
    disp(['*** MAP associations at frame ' num2str(t-min(t,Par.L)+1) ': ' num2str(assoc)]);
    
    disp(['*** Frame ' num2str(t) ' processed in ' num2str(toc) ' seconds']);
    disp('**************************************************************');
    
end

end



function [MC, BestEst, move_types] = MCMCFrame(t, L, PrevChains, PrevBest, Observs)
% Execute a frame of the fixed-lag MCMC-PF target tracker

% t - latest time frame
% L - window size
% PrevChains - MCs from the t-L to t-1 processing step, used to propose history changes
% PrevBest - The max-posterior estimate from the previous chain - used for state history <=t-L
% Observs - observations

global Par;

s = min(t,Par.S);
b = Par.BridgeLength;

% Create stores to log probabilities to prevent repeat calculations
posterior_store = -inf(Par.NumIt, Par.NumTgts);
reverse_kernel_store = -inf(Par.NumIt, Par.NumTgts);
origin_post_store = -inf(Par.NumIt, Par.NumTgts);

% Initialise the stores (the t-1 parts)
for j = 1:Par.NumTgts
    [reverse_kernel_store(1, j), ~, ~] = SampleCurrent(j, t-1, L-1, PrevBest, Observs, true);
    origin_post_store(1, j) = SingTargPosterior(j, t-1, L-1, PrevBest, Observs);
end

% Project tracks forward
for j = 1:PrevBest.N
    if t == PrevBest.tracks(j).death
        PrevBest.tracks(j) = ProjectTrack(t, PrevBest.tracks(j));
    end
end

% Initialise the stores (the t parts)
for j = 1:Par.NumTgts
    posterior_store(1, j) = SingTargPosterior(j, t, L, PrevBest, Observs);
end

% MC = Chain(Par.NumIt, PrevBest);
MC = struct( 'particles', [], 'posteriors', [], 'accept', [], 'd_accept', [] );
MC.particles = cell(Par.NumIt, 1);
MC.posteriors = -inf(Par.NumIt, Par.NumTgts);
MC.particles{1} = PrevBest;

% Create some diagnostics arrays to log move types and acceptances
accept = zeros(4,1);
d_accept = zeros(L,1);
move_types = zeros(Par.NumIt,1);

% Loop through iterations
for ii = 2:Par.NumIt
    
    % Copy the previous estimates
    Old = MC.particles{ii-1};
    New = Old;
    
    % Copy the probability stores
    MC.posteriors(ii,:) = MC.posteriors(ii-1,:);
    posterior_store(ii,:) = posterior_store(ii-1,:);
    reverse_kernel_store(ii,:) = reverse_kernel_store(ii-1,:);
    origin_post_store(ii,:) = origin_post_store(ii-1,:);
    
    % Restart chain if required
    if (mod(ii, Par.Restart)==1)
        
        k = max(1,t-1);
        weights = exp(sum(PrevChains{k}.posteriors, 2));
        weights = weights / sum(weights);
        new_part = randsample(size(PrevChains{k}.particles, 1), 1, true, weights);
        Old = PrevChains{k}.particles{new_part};
        
        % Project tracks forward
        for j = 1:Old.N
            if t == Old.tracks(j).death
                Old.tracks(j) = ProjectTrack(t, Old.tracks(j));
            end
        end
        
        New = Old;
        
        for j = 1:Par.NumTgts
            MC.posteriors(ii-1, j) = -inf;
            MC.posteriors(ii, j) = -inf;
            posterior_store(ii-1, j) = SingTargPosterior(j, t, L, New, Observs);
            reverse_kernel_store(ii-1, j) = SampleCurrent(j, t-1, L-1, New, Observs, true);
            origin_post_store(ii-1, j) = SingTargPosterior(j, t-1, L-1, New, Observs);
        end
        
    end
    
    % Choose target
%     j = unidrnd(New.N);
    j = mod(ii-2, New.N)+1;
    
    % Randomly select move type
    if t > Par.L
        type_weights = [2 1 1];
%         type_weights = [1 0 0];
    else
        type_weights = [1 0 0];
    end
    type = randsample(1:length(type_weights), 1, true, type_weights);
    move_types(ii) = type;
    
    % Switch on move type
    switch type
        
        case 1 % Single target, fixed-history
            
            % Choose proposal start-point
            d = unidrnd(L);
%             d = L;
            
            % Sample proposal
            [new_ppsl, assoc, state, mean, var] = SampleCurrent(j, t, d, New, Observs, false);
            New.tracks(j).state(t-d+1 -New.tracks(j).birth+1:t -New.tracks(j).birth+1) = state;
            New.tracks(j).covar(t-d+1 -New.tracks(j).birth+1:t -New.tracks(j).birth+1) = var;
            New.tracks(j).assoc(t-d+1 -New.tracks(j).birth+1:t -New.tracks(j).birth+1) = assoc;
            [old_ppsl, ~, ~] = SampleCurrent(j, t, d, Old, Observs, true);
            
            
        case 2 % Single target, history and window - assumes independence for frames <= t-L
            
%             sn = s;
            sn = unidrnd(s);
            k = t-sn;
            
            % Propose a new origin and copy it
            new_part = unidrnd(size(PrevChains{k}.particles, 1));
            New.tracks(j) = PrevChains{k}.particles{new_part}.tracks(j);
            NewOrigin = New;
            
            New.origin(j) = new_part;
            New.origin_time(j) = t-sn;
            
            % Extend track up to time t with blanks
            New.tracks(j).state = [New.tracks(j).state; repmat({zeros(4,1)}, sn, 1)];
            New.tracks(j).smooth = [New.tracks(j).smooth; repmat({zeros(4,1)}, sn, 1)];
            New.tracks(j).covar = [New.tracks(j).covar; repmat({zeros(4)}, sn, 1)];
            New.tracks(j).assoc = [New.tracks(j).assoc; zeros(sn, 1)];
            New.tracks(j).death = t+1;
            New.tracks(j).num = New.tracks(j).death - New.tracks(j).birth;
            
            % Propose associations
            [new_ppsl, assoc, state, mean, var] = SampleCurrent(j, t, L, New, Observs, false);
            New.tracks(j).state(t-L+1 -New.tracks(j).birth+1:t -New.tracks(j).birth+1) = state;
            New.tracks(j).covar(t-L+1 -New.tracks(j).birth+1:t -New.tracks(j).birth+1) = var;
            New.tracks(j).assoc(t-L+1 -New.tracks(j).birth+1:t -New.tracks(j).birth+1) = assoc;
            old_ppsl = SampleCurrent(j, t, L, Old, Observs, true);
            
            
        case 3 % Single target, history and bridging region - assumes independence for frames <= t-L
            
%             sn = s;
            sn = unidrnd(s);
            k = t-sn;
            
            % Propose a new origin
            new_part = unidrnd(size(PrevChains{k}.particles, 1));
            NewOrigin = New;
            NewOrigin.tracks(j) = PrevChains{k}.particles{new_part}.tracks(j);
            
            % Copy history
            ori_cut_pt = t-L - NewOrigin.tracks(j).birth + 1;
            des_cut_pt = t-L - New.tracks(j).birth + 1;
            New.tracks(j).state = [NewOrigin.tracks(j).state(1:ori_cut_pt); New.tracks(j).state(des_cut_pt+1:end)];
            New.tracks(j).smooth = [NewOrigin.tracks(j).smooth(1:ori_cut_pt); New.tracks(j).smooth(des_cut_pt+1:end)];
            New.tracks(j).assoc = [NewOrigin.tracks(j).assoc(1:ori_cut_pt); New.tracks(j).assoc(des_cut_pt+1:end)];
            New.tracks(j).birth = NewOrigin.tracks(j).birth;
            New.tracks(j).num = New.tracks(j).death - New.tracks(j).birth;
            
            New.origin(j) = new_part;
            New.origin_time(j) = t-sn;
            
            % Propose associations
            [new_ppsl, assoc, state, mean, var] = SampleCurrent(j, t-L+b, b, New, Observs, false);
            New.tracks(j).state(t-L+1 -New.tracks(j).birth+1:t-L+b -New.tracks(j).birth+1) = state;
            New.tracks(j).covar(t-L+1 -New.tracks(j).birth+1:t-L+b -New.tracks(j).birth+1) = var;
            New.tracks(j).assoc(t-L+1 -New.tracks(j).birth+1:t-L+b -New.tracks(j).birth+1) = assoc;
            old_ppsl = SampleCurrent(j, t-L+b, b, Old, Observs, true);
            
    end
    
    % Find origin posteriors and reverse kernels
    if ~(type==1)
        new_origin_post = SingTargPosterior(j, t-sn, L-sn, NewOrigin, Observs);
        if (sn < L)
            new_reverse_kernel = SampleCurrent(j, t-sn, L-sn, NewOrigin, Observs, true);
        else
            new_reverse_kernel = 0;
        end
    else
        new_origin_post = origin_post_store(ii-1,j);
        new_reverse_kernel = reverse_kernel_store(ii-1,j);
    end
    
    old_origin_post = origin_post_store(ii-1,j);
    old_reverse_kernel = reverse_kernel_store(ii-1,j);
    
    % Calculate posteriors
    new_post = SingTargPosterior(j, t, L, New, Observs);
    old_post = posterior_store(ii-1, j);
    
    % Test for acceptance
    ap = (new_post - old_post) ...
        + (old_origin_post - new_origin_post) ...
        + (old_ppsl - new_ppsl) ...
        + (new_reverse_kernel - old_reverse_kernel);
    
    if old_post==-inf, ap = inf; end
    if new_post==-inf, ap = -inf; end
    if (old_ppsl==-inf)&&(ii==2), ap = inf; end
    
    if isnan(ap)
        ap = -inf;
    end

    if log(rand) < ap
        
        % If RB, smooth state
        if Par.FLAG_RB
            last = min(t, New.tracks(j).death - 1);
            first = max(1, New.tracks(j).birth+1);
            num = last - first + 1;
            Obs = ListAssocObservs(last, num, New.tracks(j), Observs);
            init_state = New.tracks(j).state{1 -New.tracks(j).birth+1};
            init_var = Par.KFInitVar*eye(4);
            [Mean, ~] = KalmanSmoother(Obs, init_state, init_var);
            New.tracks(j).smooth(first -New.tracks(j).birth+1:last -New.tracks(j).birth+1) = Mean;
        end
     
        MC.particles{ii} = New;
        MC.posteriors(ii,j) = new_post;
        accept(type) = accept(type) + 1;
        
        posterior_store(ii,j) = new_post;
        reverse_kernel_store(ii,j) = new_reverse_kernel;
        origin_post_store(ii,j) = new_origin_post;
        
        if type==1
            d_accept(d) = d_accept(d) + 1;
        end
   
    else
        MC.particles{ii} = Old;
        MC.posteriors(ii,j) = old_post;
        
        posterior_store(ii,j) = old_post;
        %         ppsl_store(ii,j) = old_ppsl;
        reverse_kernel_store(ii,j) = old_reverse_kernel;
        origin_post_store(ii,j) = old_origin_post;
        
    end
    
end

% Pick the best particle
total_post = sum(MC.posteriors, 2);
best_ind = find(total_post==max(total_post), 1);
BestEst = MC.particles{best_ind};
BestEst.origin(:) = best_ind;
BestEst.origin_time(:) = t;

% for j = 1:Par.NumTgts
%     best_ind = find(MC.posteriors(:,j)==max(MC.posteriors(:,j)), 1);
%     BestEst.tracks(j) = MC.particles{best_ind}.tracks(j);
%     BestEst.origin(j) = best_ind;
% end

MC.accept = accept;
MC.d_accept = d_accept;

disp(['*** Accepted ' num2str(accept(1)) ' of ' num2str(sum(move_types==1)) ' fixed-history single target moves in this frame, which by d values: ' num2str(d_accept')]);
disp(['*** Accepted ' num2str(accept(2)) ' of ' num2str(sum(move_types==2)) ' full single target moves in this frame']);
% disp(['*** Accepted ' num2str(accept(3)) ' of ' num2str(sum(move_types==3)) ' full with preserved window-history single target moves in this frame']);
disp(['*** Accepted ' num2str(accept(3)) ' of ' num2str(sum(move_types==3)) ' bridging-history single target moves in this frame']);

end

