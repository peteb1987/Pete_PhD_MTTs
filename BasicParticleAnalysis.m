function [RMSE, LostTracks] = BasicParticleAnalysis(TrueTracks, Results)
%BASICANALYSIS Calculate basic performance stats for a particle-based
% tracker output.

global Par;

SE = zeros(Par.T, Par.NumTgts);
LostTracks = 0;

for t = Par.L:Par.T
    
    for j = 1:Par.NumTgts
        
%         % Find modal state
%         assoc = cellfun(@(x) x.tracks(j).assoc(t-Par.L+1 -x.tracks(j).birth+1), Results{t}.particles);
%         mode_ass = mode(assoc);
%         
%         % Find corresponding states
%         state = cellfun(@(x) x.tracks(j).state(t-Par.L+1 -x.tracks(j).birth+1), Results{t}.particles);
%         state(assoc~=mode_ass)=[];
%         
%         % Average them
%         est_state = mean(cell2mat(state'),2);
        
%         % Get MAP state
%         MAP_idx = find(Results{t}.posteriors==max(Results{t}.posteriors), 1);
%         est_state = Results{t}.particles{MAP_idx}.tracks(j).state{t-Par.L+1 -Results{t}.particles{MAP_idx}.tracks(j).birth+1};

        % MMSE estimate
        state = cellfun(@(x) x.tracks(j).state(t-Par.L+1 -x.tracks(j).birth+1), Results{t}.particles);
        est_state = mean(cell2mat(state'),2);

        % Get true state
        true_state = TrueTracks{j}.state{t-Par.L+1 -TrueTracks{j}.birth+1};
        
        % Squared Error
        SE(t, j) = sum((true_state(1:2)-est_state(1:2)).^2);
        
    end
    
end

SE(1:Par.L-1, :) = [];
RMSE = sqrt(mean(SE(:)));

end

