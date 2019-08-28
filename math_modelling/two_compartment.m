function [sumRes, S] = two_compartment(x, true_S, TEs)
% Single voxel case

S0 = (x(1))^2; % scalar
T2_1 = 1000*(sin(x(2)))^2; % scalar
T2_2 = 1000*(sin(x(3)))^2; % scalar
m1 = (sin(x(4)))^2; % scalar
m2 = 1-m1; % scalar

S = S0* ( m1*exp(-TEs/ T2_1) +  m2*exp(-TEs/ T2_2) ); % (1, num_TEs)

% Compute the sum of square differences
sumRes = sum((true_S - S).^2);