function [nodes, weights] = gauss_laguerre_quadrature(level, pdf_settings)
%   GAUSS_LAGUERRE_QUADRATURE determines the nodes and weights for the
%   Gauss-LAGUERRE rule
%
%   nodes = GAUSS_LAGUERRE_QUADRATURE(level) computes the nodes for the
%   gauss-hermite rule corresponding to quadrature order level
%   [nodes, weights] = GAUSS_LAGUERRE_QUADRATURE(level) computes the nodes
%   and corresponding weights for the gauss-laguerre rule corresponding to
%   quadrature order level
%   
%   INPUT: 
%           level: Quadrature order for which to compute the quadrature
%           rule. The number of points returned is equal to level [Integer]
%           alpha: Alpha setting for Laguerre polynomial [float]
%
%   OUTPUT: 
%           nodes: Nodes for the gauss-laguerre rule [level X 1 vector]
%           weights: Weights for the gauss-laguerre rule [level X 1 vector]
%
%   See also LAGUERRE

    
    % Constructing the companion matrix based on the level, just need to
    % create the diagonals based on the recurrence relation
    i = 1:level;
    a = 2.*i - 1;
    b = abs(i);
    % Cut off last element of b, don't need it
    b = b(1:level-1);
    
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
    
    nodes = nodes.*pdf_settings.scale + pdf_settings.a;

    % The corresponding (normalized) weights are just the first element of the eigen
    % vector squared.
    weights = transpose(eigen_vec(1,index).^2);
end