function [sumRes, S] = singleCompartment(x, true_S, TEs)
% Single voxel case

S0 = (x(1))^2; % scalar
T2 = (x(2))^2; % scalar

if T2>200
  T2=200;
end

S = S0*exp(-TEs/ T2);

% Compute the sum of square differences
sumRes = sum((true_S - S).^2);