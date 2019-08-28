function logL = logLikelihood(data, model, parnames, parvals)

% get the data from the data cell array
TEs = data{1};
A = data{2};
if length(data)==3
    noise_sd = data{3};
else
    noise_sd = 0;
end

% evaluate the model
S = feval(model, TEs, parnames, parvals);
% if the model returns a NaN then set the likelihood to be zero (e.g. loglikelihood to be -inf)
if isnan(S)
    logL = -inf;
    return;
end

% calculate the log likelihood
if noise_sd ~= 0
    logL = log(sum((A-S').^2)/(2*noise_sd^2));
else
    logL = mean(log(sum((A-S).^2)));
end

if isnan(logL)
    error('Error: log likelihood is NaN!');
end

end

