function [Resampled] = ChainSystematicResample(Chain, weight )
%SYSTEMATICRESAMPLE Resample a single target systematically

global Par;

Np = Par.NumIt;

% Change the weights to linear
weight = exp(weight);

% Copy distribution
Resampled = Chain;

% Create empty vector for offspring count
N = zeros(Np, 1);

% Generate random index array
u = zeros(Np, 1);
u(1) = rand/Np;
for ii = 2:Np
    u(ii) = u(1) + (ii-1)/Np;
end    

% Generate cumulative weight array
w_sum = cumsum(weight);

% Enumerate offspring
for ii = 1:Np
    if ii>1
        N(ii) = sum((u < w_sum(ii))&(u > w_sum(ii-1)));
    else
        N(ii) = sum(u < w_sum(ii));
    end
end

assert(sum(N)==Np, 'Wrong number of resampled children');

% Construct new, resampled array
jj = 1;
for ii = 1:Np
    for k = 1:N(ii)
        Resampled.particles{jj} = Chain.particles{ii};
        jj = jj + 1;
    end
end

Resampled.weights = log(ones(Np, 1)/Np);

end