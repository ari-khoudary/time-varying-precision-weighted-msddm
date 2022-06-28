function variance = betaVar(alpha, beta)

% computes the variance of a beta distribution

variance = (alpha*beta) / ((alpha+beta)^2) * (alpha + beta + 1);

end