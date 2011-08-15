function [RMSE_MMSE, RMSE_MAP, prop_ass, num_lost] = BasicParticleAnalysis(TrueTracks, Results)
%BASICANALYSIS Calculate basic performance stats for a particle-based
% tracker output.

global Par;

SE_MMSE = zeros(Par.T, Par.NumTgts);
SE_MAP = zeros(Par.T, Par.NumTgts);
LostTracks = false(1, Par.NumTgts);
CorrectAssociations = zeros(Par.T, Par.NumTgts);

for t = Par.AnalysisLag:Par.T
    
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
        
        % Get MAP state
        MAP_idx = find(sum(Results{t}.posteriors,2)==max(sum(Results{t}.posteriors,2)), 1);
        if ~Par.FLAG_RB
            MAP_state = Results{t}.particles{MAP_idx}.tracks(j).state{t-Par.AnalysisLag+1 -Results{t}.particles{MAP_idx}.tracks(j).birth+1};
        else
            MAP_state = Results{t}.particles{MAP_idx}.tracks(j).smooth{t-Par.AnalysisLag+1 -Results{t}.particles{MAP_idx}.tracks(j).birth+1};
        end

        % MMSE estimate
        if ~Par.FLAG_RB
            state = cellfun(@(x) x.tracks(j).state(t-Par.AnalysisLag+1 -x.tracks(j).birth+1), Results{t}.particles);
        else
            state = cellfun(@(x) x.tracks(j).smooth(t-Par.AnalysisLag+1 -x.tracks(j).birth+1), Results{t}.particles);
        end
        MMSE_state = mean(cell2mat(state'),2);

        % Get true state
        true_state = TrueTracks{j}.state{t-Par.AnalysisLag+1 -TrueTracks{j}.birth+1};
        
        % Squared Error
        SE_MMSE(t, j) = sum((true_state(1:2)-MMSE_state(1:2)).^2);
        SE_MAP(t, j) = sum((true_state(1:2)-MAP_state(1:2)).^2);
        
        % Count correct associations
        assoc = cellfun(@(x) x.tracks(j).assoc(t-Par.AnalysisLag+1 -x.tracks(j).birth+1), Results{t}.particles);
        true_assoc = TrueTracks{j}.assoc(t-Par.AnalysisLag+1 -TrueTracks{j}.birth+1);
        CorrectAssociations(t,j) = sum(assoc==true_assoc);
        
        % Find lost tracks
        if CorrectAssociations(t,j) == 0
            tracking = false;
        else
            tracking = true;
        end
        for tt = t-Par.AnalysisLag+2:t
            assoc = cellfun(@(x) x.tracks(j).assoc(tt -x.tracks(j).birth+1), Results{t}.particles);
            true_assoc = TrueTracks{j}.assoc(tt -TrueTracks{j}.birth+1);
            if (true_assoc~=0)&&(sum(assoc==true_assoc)>0)
                tracking = true;
                break;
            end
        end
        
        if ~tracking
            LostTracks(j) = true;
        end
        
    end
    
end

prop_ass = sum(CorrectAssociations(:)) / ( Par.NumTgts*(Par.T-Par.AnalysisLag+1)*length(Results{end}.particles) );
num_lost = sum(LostTracks);

SE_MMSE(1:Par.AnalysisLag-1, :) = [];
SE_MAP(1:Par.AnalysisLag-1, :) = [];
SE_MMSE(:,LostTracks) = [];
SE_MAP(:,LostTracks) = [];
RMSE_MMSE = sqrt(mean(SE_MMSE(:)));
RMSE_MAP = sqrt(mean(SE_MAP(:)));

end

