function [ Observs, hits ] = GenerateObservations( TrueTracks )
%GENERATEOBSERVATIONS Generate a set of observations from a set of known tracks

global Par;

% Simple Observation Generator - Observe target with probability PDetect and
% Gaussian observation noise ObsNoiseVar. Poisson clutter.

% Make some local variables for convenience
T = Par.T;               % Number of frames

% Create a structure to store observations
obs = struct('r', [], 'N', []);             % r is an array of observations (one per row).
                                            % N is the number of observations
Observs = repmat(obs, Par.T, 1);

% Create a cell array to store inidices of target observations
hits = zeros(T, Par.NumTgts);

% Loop over time
for t = 1:T
    
    % Draw number of clutter obs from a poisson dist
    num_clut = poissrnd(Par.ExpClutObs);
    
    % Count the number of targets present in this frame
    num_tgts = 0;
    for j = 1:Par.NumTgts
        if (TrueTracks{j}.birth <= t)&&(TrueTracks{j}.death > t)
            num_tgts = num_tgts + 1;
        end
    end
    
    Observs(t).N = num_tgts+num_clut;
    
    % Initialise cell array for time instant with max required size
    Observs(t).r = zeros(Observs(t).N, 2);
    
    % Initialise observation counter
    i = 1;
    
    % Generate target observations
    for j = 1:Par.NumTgts
        if (TrueTracks{j}.birth <= t)&&(TrueTracks{j}.death > t)
            state = TrueTracks{j}.state{t-TrueTracks{j}.birth+1};
            
            hits(t, j) = i;
            
            if Par.FLAG_ObsMod == 0
                % Cartesian observations with gaussian noise
                Observs(t).r(i, :) = mvnrnd(state(1:2)', Par.ObsNoiseVar*ones(1,2));
                
            elseif Par.FLAG_ObsMod == 1
                % Polar observations with gaussian noise
                [bng, rng] = cart2pol(state(1), state(2));
                Observs(t).r(i, :) = mvnrnd([bng, rng], Par.R);
            end
            
            % Remove missed detections
            if rand < (1-Par.PDetect)
                hits(t, j) = 0;
                Observs(t).r(i, :) = [];
                Observs(t).N = Observs(t).N - 1;
                num_tgts = num_tgts - 1;
                i = i - 1;
            end
            
            i = i + 1;
        end
    end
    
    % Generate clutter observations
    for i = num_tgts+1:Observs(t).N
        
        if Par.FLAG_ObsMod == 0
            % Cartesian Poisson clutter
            Observs(t).r(i, 1) = unifrnd(-Par.Xmax, Par.Xmax);
            Observs(t).r(i, 2) = unifrnd(-Par.Xmax, Par.Xmax);

        elseif Par.FLAG_ObsMod == 1
            % Polar Poisson clutter
            Observs(t).r(i, 1) = unifrnd(-pi, pi);
            Observs(t).r(i, 2) = unifrnd(0, Par.Xmax);
            
        end
    end
    
end


end

