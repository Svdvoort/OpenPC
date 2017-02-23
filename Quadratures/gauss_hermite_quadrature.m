function [nodes, weights] = gauss_hermite_quadrature(level, pdf_settings)
%   GAUSS_HERMITE_QUADRATURE determines the nodes and weights for the
%   Gauss-Hermite rule
%
%   nodes = GAUSS_HERMITE_QUADRATURE(level) computes the nodes for the
%   gauss-hermite rule corresponding to quadrature order level
%   [nodes, weights] = GAUSS_HERMITE_QUADRATURE(level) computes the nodes
%   and corresponding weights for the gauss-hermite rule corresponding to
%   quadrature order level
%   
%   INPUT: 
%           level: Quadrature order for which to compute the quadrature
%           rule. The number of points returned is equal to level [Integer]
%
%   OUTPUT: 
%           nodes: Nodes for the gauss-hermite rule [level X 1 vector]
%           weights: Weights for the gauss-hermite rule [level X 1 vector]
%
%   In this function probabilist' Hermite polynomials have been used, which
%   will result in different nodes and weights then when physicists'
%   Hermite polynomials are used. The probablist version is used as these
%   are monic polynomials. 
%
%   See also HERMITE

    % Constructing the companion matrix based on the level, just need to
    % create the diagonals based on the recurrence relation
    b = sqrt(1:(level-1));
    companion_mat = diag(b, -1) + diag(b, 1);
    
    if verLessThan('matlab', '8.3')
        % Exception for matlab release before 2014
        [eigen_vec, eigen_val] = eig(companion_mat);
        eigen_val = diag(eigen_val);
    else      
        [eigen_vec, eigen_val] = eig(companion_mat, 'vector');
    end
    
    % Nodes are eigenvalues of the companion matrix
    [nodes, index] = sort(eigen_val);
    
    % If there is an uneven number of points, the middle point should
    % actually be 0 by definition, but small deviations occur  because of
    % computational limitations
    if mod(level, 2) ~= 0
        nodes((level+1)/2) = 0;
    end
    
    % Scale the nodes according to the pdf settings
    nodes = nodes*pdf_settings.std + pdf_settings.mean;
    
    % The corresponding (normalized) weights are just the first element of the eigen
    % vector squared.
     weights = transpose(eigen_vec(1, index).^2);
end