function [nodes, weights] = arbitrary_quadrature(level, pdf_settings)

    % Constructing the companion matrix based on the level, just need to
    % create the diagonals based on the recurrence relation
    
    % First element is excluded as this is beta_0 and gamma_0, which should
    % not be in matrix
    
    beta = pdf_settings.beta;
    gamma = pdf_settings.gamma;
    a = -beta(1:level);
    b = sqrt(gamma(2:level));

    companion_mat = diag(b, -1) + diag(a, 0) + diag(b, 1);
        
    if verLessThan('matlab', '8.3')
        % Exception for matlab release before 2014
        [eigen_vec, eigen_val] = eig(companion_mat);
        eigen_val = diag(eigen_val);
    else      
        [eigen_vec, eigen_val] = eig(companion_mat, 'vector');
    end    
    
    % Nodes are eigenvalues of the companion matrix
    [nodes, index] = sort(eigen_val);
        
    % The corresponding (normalized) weights are just the first element of the eigen
    % vector squared.
     weights = transpose(eigen_vec(1, index).^2);
end