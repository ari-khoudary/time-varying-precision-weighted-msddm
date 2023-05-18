function variance = betaVar(alpha, beta)

% computes the variance of a beta distribution
variance = zeros(length(alpha), 1);
variance(:, 1) = (alpha(:, 1).*beta(:, 1)) ./ (((alpha(:, 1)+beta(:,1)).^2) .* (alpha(:,1) + beta(:,1) + 1));

end