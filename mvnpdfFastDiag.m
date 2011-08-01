function [ pdf ] = mvnpdfFastDiag( X, Mu, Var )
%MVNPDFFASTSYMM This is an alternative to mvnpdf when the covariance matrix
%is diagonal. Var is a column vector.

A = 0.5*(-sum(((X-Mu).^2)./Var)-log(prod(Var))-length(X)*log(2*pi));

pdf = exp(A);

end