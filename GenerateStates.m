function [ TrueTracks, InitStates ] = GenerateStates( Accelerations )
%GENERATESTATES Generate target states. Parameters, including which model
%to use, are in global Par.

global Par;
global Templates;

if nargin < 1
    Accelerations = repmat({zeros(Par.T, 2)}, Par.NumTgts, 1);
end

% Number of targets
N = Par.NumTgts;

% Number of time steps
T = Par.T;

% Initialise a cell array for the states
TrueTracks = cell(Par.NumTgts, 1);
InitStates = cell(Par.NumTgts, 1);

% Set default parameters
for j = 1:N
    % Create track structure
    TrueTracks{j} = Templates.Track;
    
    % Add basic properties
    TrueTracks{j}.birth = 0;
    TrueTracks{j}.death = T + 1;
    TrueTracks{j}.num = TrueTracks{j}.death - TrueTracks{j}.birth;
    
    % Initialise state and assoc arrays
    TrueTracks{j}.state = repmat({zeros(4,1)}, TrueTracks{j}.num, 1);
    TrueTracks{j}.assoc = zeros(TrueTracks{j}.num, 1);
    
    % Randomly generate initial state
    if Par.FLAG_ObsMod == 0
        TrueTracks{j}.state{1}(1:2) = unifrnd(-Par.MaxInitStateDist*Par.Xmax, Par.MaxInitStateDist*Par.Xmax, [2,1]);
    elseif Par.FLAG_ObsMod == 1
         rng = unifrnd(Par.MinInitStateRadius*Par.Xmax, Par.MaxInitStateRadius*Par.Xmax);
         bng = unifrnd(0, 2*pi);
         TrueTracks{j}.state{1}(1) = rng*sin(bng);
         TrueTracks{j}.state{1}(2) = rng*cos(bng);
    end
    TrueTracks{j}.state{1}(3:4) = unifrnd(-Par.Vmax, Par.Vmax, [2,1]);
    
    if Par.FLAG_SetInitStates
        TrueTracks{j}.state{1} = Par.InitStates{j};
    end
    
    if TrueTracks{j}.birth == 0;
        InitStates{j} = TrueTracks{j}.state{1};
    end
    
    % Loop through frames
    for k = 2:TrueTracks{j}.num

        % Calculate expected state
        exp_state = Par.A * TrueTracks{j}.state{k-1} + Par.B * Accelerations{j}(k-1, :)';
        
        % Sample state from Gaussian
        TrueTracks{j}.state{k} = mvnrnd(exp_state', Par.Q)';

%         % Kill if outside scene
%         if (Par.FLAG_ObsMod==0) && any(abs( TrueTracks{j}.state{k}(1:2))>Par.Xmax) || ...
%                 (Par.FLAG_ObsMod==1) && any(sqrt(sum( TrueTracks{j}.state{k}(1:2).^2))>Par.Xmax)
%             TrueTracks{j}.state(k:end) = [];
%             TrueTracks{j}.assoc(k:end) = [];
%             TrueTracks{j}.death = TrueTracks{j}.birth + k - 1;
%             TrueTracks{j}.num = TrueTracks{j}.death - TrueTracks{j}.birth;
%             break;
%         end
        
        % Limit velocity
        TrueTracks{j}.state{k}(3:4) = min( max( TrueTracks{j}.state{k}(3:4), -Par.Vmax), Par.Vmax);
        
    end
    
end

end

