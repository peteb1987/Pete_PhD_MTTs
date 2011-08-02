function [ obs ] = ListAssocObservs( t, L, Track, Observs )
%LISTASSOCOBSERVS Draw up a list of observations associated with Track
% between t-L+1 and t.

obs = cell(L, 1);

% Loop through window
for k = 1:L
    
    % Get the association index
    ass = Track.assoc(t-L+k -Track.birth+1);
    
    % Get observation if one is associated
    if ass ~= 0
        obs{k} = Observs(t-L+k).r(ass, :)';
    end
    
end

end

