function [sumRes, S] = three_compartment(x, true_S, TEs)
% Single voxel case

S0 = (x(1))^2; % scalar
m1 = sin(x(2))^2;
m2 = (1-m1)*sin(x(3))^2;
m3 = 1-m1-m2;

S = S0 * ( m1*exp(-TEs/20) + m2*exp(-TEs/80) + m3*exp(-TEs/2000)); % (1, num_TEs)

% Compute the sum of square differences
sumRes = sum((true_S - S).^2);
end

