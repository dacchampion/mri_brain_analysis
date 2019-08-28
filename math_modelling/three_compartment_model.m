function S = three_compartment_model(TEs, parnames, parvals)
    lpn = length(parnames);
    lpv = length(parvals);
    if lpn ~= lpv
        error('Error: Number of parameters and their values should be the same size!');
    end
    
    % extract parameter values
    S0 = parvals{find(strcmp(parnames, 'S0'))};
    m1 = parvals{find(strcmp(parnames, 'm1'))};
    m2 = parvals{find(strcmp(parnames, 'm2'))}; 
    m3 = parvals{find(strcmp(parnames, 'm3'))};
    S = S0 * (m1*exp(-TEs/20) + m2*exp(-TEs/80) + m3*exp(-TEs/2000)); 
end

