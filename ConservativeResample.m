function [Resampled] = ConservativeResample(j, PartSet, weight )
%CONSERVATIVERESAMPLE Resample a single target conservatively

global Par;

Np = Par.NumPart;

% Change the weights to linear
weight = exp(weight);

% Copy distribution
Resampled = PartSet;

% Create empty vector for offspring count
Nchild = zeros(Np, 1);

% Enumerate offspring
for ii = 1:Np
    Nchild(ii) = max(1,floor(  1 *  weight(ii)*Np));
end
Nct = sum(Nchild);

% Generate an array of selection indices and reweight
parent = zeros(Nct,1);
interm_weight = zeros(Nct,1);
jj=1;
for ii = 1:Np
    w = weight(ii)/Nchild(ii);
    parent(jj:jj+Nchild(ii)-1) = ii;
    interm_weight(jj:jj+Nchild(ii)-1) = w;
    jj = jj+Nchild(ii);
end

% Get rid of zero weights
parent(interm_weight==0)=[];
interm_weight(interm_weight==0)=[];

% Get rid of low weights
thresh = max(interm_weight) / exp(Par.ResampleLowWeightThresh);
parent(interm_weight<thresh)=[];
interm_weight(interm_weight<thresh)=[];
Nct = length(parent);

% Systematically downsample
u = ceil( Nct * cumsum( (1/Np)*ones(Np,1) ) - unifrnd(0, 1/Np) );
samp = parent(u);
new_weight = interm_weight(u);

new_weight = new_weight/sum(new_weight);
new_weight = log(new_weight);

% Construct new, resampled array
for ii = 1:Np
    Resampled.particles{ii}.tracks(j) = PartSet.particles{samp(ii)}.tracks(j);
end
Resampled.weights{j} = new_weight;

end