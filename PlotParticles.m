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

global Par;

% Loop through targets
for j = 1:Part.N
    
    if Part.tracks(j).num > 0
        
        % Choose a colour
        col = [rand, rand, 0];
        
        % Collate state
        if ~Par.FLAG_RB
            state = cell2mat(Part.tracks(j).state');
        else
            state = cell2mat(Part.tracks(j).smooth');
        end
        x = state(1, :);
        y = state(2, :);
        
        % Plot track
        plot(x, y, '-', 'color', col);
        
        if (Par.FLAG_AlgType == 2)||(Par.FLAG_AlgType == 3)
            state = Part.tracks(j).state;
            covar = Part.tracks(j).covar;
            for t = 1:length(state)
                plot_gaussian_ellipsoid(state{t}(1:2), covar{t}(1:2,1:2));
            end
        end
        
    end
    
end

end

