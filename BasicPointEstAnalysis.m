function [RMSE, num_lost] = BasicPointEstAnalysis(TrueTracks, Results)
%BASICANALYSIS Calculate basic performance stats for a particle-based
% tracker output.

global Par;

SE = zeros(Par.T, Par.NumTgts);
LostTracks = false(1, Par.NumTgts);


for t = Par.L:Par.T
    
    for j = 1:Par.NumTgts
        
        % MMSE estimate
        state = Results{t}.tracks(j).state{t-Par.L+1 -Results{t}.tracks(j).birth+1};

        % Get true state
        true_state = TrueTracks{j}.state{t-Par.L+1 -TrueTracks{j}.birth+1};
        
        % Squared Error
        SE(t, j) = sum((true_state(1:2)-state(1:2)).^2);
        
        % Find lost tracks
        tracking = false;
        for tt = t-min(t,5)+1:t
            if (SE(t,j)<100)
                tracking = true;
                break;
            end
        end
        
        if ~tracking
            LostTracks(j) = true;
        end
        
    end
    
end

num_lost = sum(LostTracks);

SE(1:Par.L-1, :) = [];
SE(:,LostTracks) = [];
RMSE = sqrt(mean(SE(:)));

end

