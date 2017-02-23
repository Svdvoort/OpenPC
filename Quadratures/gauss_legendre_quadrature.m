function [nodes, weights] = gauss_legendre_quadrature(level, pdf_settings)
%   GAUSS_LEGENDRE_QUADRATURE determines the nodes and weights for the
%   Gauss-Legendre rule
%
%   nodes = GAUSS_LEGENDRE_QUADRATURE(level) computes the nodes for the
%   gauss-legendre rule corresponding to quadrature order level
%   [nodes, weights] = GAUSS_LEGENDRE_QUADRATURE(level) computes the nodes
%   and corresponding weights for the gauss-legendre rule corresponding to
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
%   See also LEGENDRE

    % Constructing the companion matrix based on the level, just need to
    % create the diagonals based on the recurrence relation
    i = 1:level-1;
    b = sqrt(i.^2./(4.*i.^2 -1));
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

    pdf_range = (pdf_settings.b - pdf_settings.a)/2;
    pdf_middle = (pdf_settings.b + pdf_settings.a)/2;
    
    nodes = nodes * pdf_range + pdf_middle;
    % The corresponding (normalized) weights are just the first element of the eigen
    % vector squared.
    weights = transpose(eigen_vec(1, index).^2);
end