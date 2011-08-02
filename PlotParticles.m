function f = PlotParticles( Particles, f )
%PLOTPARTICLES Plots all the tracks contained in Particles

global Par;

if nargin == 1
    f = figure; hold on
    xlim([-Par.Xmax Par.Xmax]), ylim([-Par.Xmax Par.Xmax])
else
    figure(f);
end

cellfun(@PlotOnePart, Particles);

plot(0, 0, 'xk');

end

function PlotOnePart(Part)

% Loop through targets
for j = 1:Part.N
    
    if Part.tracks(j).num > 0
        
        % Choose a colour
        col = [rand, rand, 0];
        
        % Collate state
        state = cell2mat(Part.tracks(j).state');
        x = state(1, :);
        y = state(2, :);
        
        % Plot track
        plot(x, y, '-', 'color', col);
        
    end
    
end

end

