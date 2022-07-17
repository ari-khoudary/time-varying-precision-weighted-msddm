function entropy = computeEntropy(probability)

if probability == 1
    probability = 0.9999999999;
elseif probability == 0
    probability = 1e-10;
end

entropy = -probability*(log(probability)) - ((1-probability)*log(1-probability));

end