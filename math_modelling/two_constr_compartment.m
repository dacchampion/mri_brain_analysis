function [sumRes, S] = two_constr_compartment(x, true_S, TEs)
% Single voxel case

S0 = (x(1))^2; % scalar
m1 = (sin(x(2)))^2; % scalar
m2 = 1-m1; % scalar

S = S0* ( m1*exp(-TEs/20) +  m2*exp(-TEs/200) ); % (1, num_TEs)

% Compute the sum of square differences
sumRes = sum((true_S - S).^2);