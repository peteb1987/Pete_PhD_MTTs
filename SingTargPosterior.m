function [ post ] = SingTargPosterior(j, t, L, Set, Observs)
%SINGTARGPOSTERIOR Calculate posterior probability factor corresponging to 
% a single target. 

global Par;

% Initialise probabilities
like = zeros(L, 1);
trans = zeros(L, 1);
clut = zeros(L, 1);
assoc = zeros(L, 1);

end_time = min(t, Set.tracks(j).death-1);
start_time = max(t-L+1, Set.tracks(j).birth);

% Loop through window
for tt = start_time:end_time
    k = tt-(t-L);
    
    % See what observations have already been used
    obs_assigned = [];
    for jj = 1:Set.N
        if (j ~= jj) && ((Set.tracks(j).birth <= tt)&&(Set.tracks(j).death > tt))
            ass = Set.tracks(jj).assoc(tt -Set.tracks(jj).birth+1);
            if ass > 0
                obs_assigned = [obs_assigned, ass];
            end
        end
    end
    
    % Get states
    state = Set.tracks(j).state{tt -Set.tracks(j).birth+1};
    
    % Calculate likelihood
    if any(abs(state(3:4))>Par.Vlimit)||any(abs(state(1:2))>2*Par.Xmax)
        like(k) = -inf;
    else
        ass = Set.tracks(j).assoc(tt -Set.tracks(j).birth+1);
        if ass~=0
            if Par.FLAG_ObsMod == 0
                like(k) = like(k) + log( mvnpdfFastDiag(Observs(tt).r(ass, :), state(1:2)', diag(Par.R)') );
            elseif Par.FLAG_ObsMod == 1
                [bng, rng] = cart2pol(state(1), state(2));
                if (Observs(tt).r(ass, 1) - bng) > pi
                    bng = bng + 2*pi;
                elseif (Observs(tt).r(ass, 1) - bng) < -pi
                    bng = bng - 2*pi;
                end
                like(k) = like(k) + log( mvnpdfFastDiag(Observs(tt).r(ass, :), [bng rng], diag(Par.R)') );
            end
        end
    end
    
    % Calculate transition density
    if tt==Set.tracks(j).birth
        trans(k) = trans(k) + log(Par.UnifPosDens*Par.UnifVelDens);
    else
        prev_state = Set.tracks(j).state{tt-1 -Set.tracks(j).birth+1};
        trans(k) = trans(k) + log( (1-Par.PDeath) * mvnpdfQ(state', (Par.A * prev_state)') );
        % trans(k) = trans(k) + log( (1-Par.PDeath) * mvnpdf(state', (Par.A * prev_state)', Par.Q) );
    end
    
    % See if it's died
    if (Set.tracks(j).death == t)
        trans(k) = trans(k) + log(Par.PDeath);
    end
    
    
    % Calculate association term
    if (Set.tracks(j).birth <= tt)&&(Set.tracks(j).death > tt)
        ass = Set.tracks(j).assoc(tt -Set.tracks(j).birth+1);
        if any(ass==obs_assigned)
            assoc(k) = -inf;
            break
        elseif ass>0
            assoc(k) = assoc(k) + log(Par.PDetect);
        elseif ass==0
            assoc(k) = assoc(k) + log((1-Par.PDetect) * Par.ExpClutObs);
        else
            error('Invalid association');
        end
    else
        assoc(k) = assoc(k) + log(Par.ExpClutObs);
    end
    
    
    % Calculate clutter term
    if ass == 0
        clut(k) = log(Par.ClutDens);
    end
    
end


% if any(isinf(like))
%     disp('Zero likelihood');
% end
% if any(isinf(trans))
%     disp('Zero transition density');
% end
% if any(isinf(clut))
%     disp('Zero clutter probability');
% end
% if any(isinf(assoc))
%     disp('Zero association probability');
% end

% Combine likelihood and transition density terms
post = sum(like) + sum(trans) + sum(clut) + sum(assoc);

end