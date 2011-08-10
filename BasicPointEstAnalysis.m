function [RMSE, num_lost] = BasicPointEstAnalysis(TrueTracks, Results)
%BASICANALYSIS Calculate basic performance stats for a particle-based
% tracker output.

global Par;

SE = zeros(Par.T, Par.NumTgts);
LostTracks = false(1, Par.NumTgts);


for t = Par.AnalysisLag:Par.T
    
    for j = 1:Par.NumTgts
        
        % MMSE estimate
        state = Results{t}.tracks(j).state{t-Par.AnalysisLag+1 -Results{t}.tracks(j).birth+1};

        % Get true state
        true_state = TrueTracks{j}.state{t-Par.AnalysisLag+1 -TrueTracks{j}.birth+1};
        
        % Squared Error
        SE(t, j) = sum((true_state(1:2)-state(1:2)).^2);
        
        % Find lost tracks
        tracking = false;
        for tt = t-Par.AnalysisLag+1:t
            if (SE(t,j)<1000)
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

SE(1:Par.AnalysisLag-1, :) = [];
SE(:,LostTracks) = [];
RMSE = sqrt(mean(SE(:)));

end

