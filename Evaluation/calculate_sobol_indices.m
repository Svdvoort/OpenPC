function [sobol_indices, total_sensitivity, u_combined_terms] = calculate_sobol_indices(PCE)

    N_dims = size(PCE.basis, 2);

    
    % Compute standard deviation from coefficients
    std_dev = sum(PCE.coeffs(:, 2:end).^2 .* transpose(PCE.norm(2:end)));
    
    % Find where there are unique combined terms. Make a logical index
    % Of basis, which gives same rows for same interactions (e.g.
    % basis [1, 2, 0] is same interaction as [2, 1, 0])
    combined_terms = PCE.basis ~= 0;
    [u_combined_terms, ~, I_u_combined_terms] = unique(combined_terms, 'rows');
    N_u_combined_terms = size(u_combined_terms, 1);
    
    sobol_indices = zeros(N_u_combined_terms, 1);
    for i_u_combined_term = 1:N_u_combined_terms
        % Find index of current combined terms for which to compute sobol
        % index
        LI_combined_terms = I_u_combined_terms == i_u_combined_term;
        sobol_indices(i_u_combined_term) = PCE.coeffs(:, LI_combined_terms).^2 * PCE.norm(LI_combined_terms);
    end    
    
    sobol_indices = sobol_indices./std_dev;
    
    total_sensitivity = zeros(N_dims, 1);
    % Create matrix to contain interaction terms.
    for i_dim = 1:N_dims        
        temp = u_combined_terms(:, i_dim) == 1;
        total_sensitivity(i_dim, :) = sum(sobol_indices(temp));
    end
    

end