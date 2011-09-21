function [ sig_pts, sig_wts, Np ] = UnscentedTransform( mean, var )
%UNSCENTEDTRANSFORM Returns a set of sigma points for a mean and variance

% mean is a column vector, var is a square matrix

kappa = 1;      % Unscented transform parameter (as in Julier & Uhlmann 1997)

% Initialise
N = size(mean, 1);
Np = 2*N+1;

sig_pts = zeros(N,Np);
sig_wts = zeros(Np,1);

% Matrix square root
mat_sq_rt = chol((N+kappa)*var)';

% First point
sig_wts(1) = kappa/(N+kappa);
sig_pts(:,1) = mean;

% Positive points
for ii = 1:N
    sig_pts(:,ii+1) = mean + mat_sq_rt(:,ii);
    sig_wts(ii+1) = 1/(2*(N+kappa));
end

% Negative points
for ii = 1:N
    sig_pts(:,ii+N+1) = mean - mat_sq_rt(:,ii);
    sig_wts(ii+N+1) = 1/(2*(N+kappa));
end

end

