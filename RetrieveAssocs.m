function [ assoc, correct ] = RetrieveAssocs( fr, Particles, detections )
%RETRIEVEASSOCS Return an array of track associations for each particle

assoc = cell(Particles{1}.N, 1);
correct = [];

for j = 1:Particles{1}.N
    assoc{j} = zeros(length(Particles), fr);
    for tt = 1:fr
        assoc{j}(:,tt) = cellfun(@(x) x.tracks(j).assoc(tt+1), Particles);
    end
end

if nargin == 3
    correct = zeros(Particles{1}.N, fr);
    for j = 1:Particles{1}.N
        for t = 1:fr
            correct(j, t) = sum(assoc{j}(:,t)==detections(t,j));
        end
    end
    figure, plot(correct');
end

end

