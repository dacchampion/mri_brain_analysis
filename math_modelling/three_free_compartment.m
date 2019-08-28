function [sumRes, S] = three_free_compartment(x, true_S, TEs)
% Single voxel case

S0 = (x(1))^2; % scalar
m1 = sin(x(2))^2;
m2 = (1-m1)*sin(x(3))^2;
m3 = 1-m1-m2;
T2_1 = 20*(sin(x(4)))^2; % scalar
T2_2 = 80*(sin(x(5)))^2; % scalar
T2_3 = 2000*(sin(x(6)))^2; % scalar

S = S0 * ( m1*exp(-TEs/T2_1) + m2*exp(-TEs/T2_2) + m3*exp(-TEs/T2_3)); % (1, num_TEs)

% Compute the sum of square differences
sumRes = sum((true_S - S).^2);
end

