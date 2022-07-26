function variance = altBetaVar(mu, v)

% computes the variance of a beta distribution parameterized by mean &
% sample size

variance = (mu*(1-mu)) / (1 + v);

end