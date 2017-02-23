function quadratures = get_quadratures(grid_level, quadrature_rules, pdf_settings)
%   GET_QUADRATURE determines the nodes and weights for a certain quadrature
%
%   quadratures = GET_QUADRATURES(grid_level, quadrature_rule) gives back the
%   quadrature structure with then nodes and weights
%   
%   INPUT: 
%           grid_level: The maximum grid level for which to compute the
%           quadratures [Integer]
%           quadrature_rule: The quadrature rules to use. The currently
%           supported quadrature rules are: gauss-hermite, gauss-legendre and
%           gauss-laguerre. [N x 1 cell of strings, where N is the number of 
%           different quadratures]
%
%   OUTPUT: 
%          quadratures: The quadrature structure with a nodes and weights
%          field.

    N_quad_types = numel(quadrature_rules); 
        
    % Initialize the quadratures structure
    quadratures = struct();
    quadratures.nodes = cell(N_quad_types, grid_level);
    quadratures.weights = cell(N_quad_types, grid_level);
    quadratures.sparse_nodes = cell(N_quad_types, grid_level);
    quadratures.sparse_weights = cell(N_quad_types, grid_level);
    quadratures.quad_rules = quadrature_rules;
    
    % TODO: could still do uniqueness for everything that is not arbitrary
    % Scaling would then have to be done outside the quadrature function
    
    for i_quad_type = 1:numel(quadrature_rules)        
        switch quadrature_rules{i_quad_type}
            case 'gauss-hermite'
                quad_function_handle = @(x) gauss_hermite_quadrature(x, pdf_settings.pdf_parameters{i_quad_type});
            case 'gauss-legendre'
                quad_function_handle = @(x) gauss_legendre_quadrature(x, pdf_settings.pdf_parameters{i_quad_type});
            case 'gauss-laguerre'
                quad_function_handle = @(x) gauss_laguerre_quadrature(x, pdf_settings.pdf_parameters{i_quad_type});
            case 'gauss-jacobi'
                quad_function_handle = @(x) gauss_jacobi_quadrature(x, pdf_settings.pdf_parameters{i_quad_type});
            case 'gauss-arbitrary'
                quad_function_handle = @(x) arbitrary_quadrature(x, pdf_settings.pdf_parameters{i_quad_type});                
            otherwise
                disp(quadrature_rules{i_quad_type})
                error('get_error_scenarios:unknownQuadType','The requested quadrature type is not known');
        end
        
        [nodes, weights] = calculate_quadrature(quad_function_handle, double(grid_level));
        
        quadratures.nodes(i_quad_type, : ) = nodes;
        quadratures.weights(i_quad_type, :) = weights;
        
        % Do some pre allocation
        sparse_nodes = cell(1, grid_level);
        sparse_weights = cell(1, grid_level);
        
        sparse_nodes{1, 1} = nodes{1, 1};
        sparse_weights{1, 1} = weights{1, 1};

        for i_node = 2:size(nodes, 2)
            [sparse_nodes{1, i_node}, ~, u_node_index] = unique([nodes{1, i_node}; nodes{1, i_node-1}]);            
            sparse_weights{1, i_node} = accumarray(u_node_index, [weights{1, i_node}; -weights{1, i_node-1}]);
        end
        
        quadratures.sparse_nodes(i_quad_type, :) = sparse_nodes;
        quadratures.sparse_weights(i_quad_type, :) = sparse_weights;

    end
end


% Use an implicit function, as otherwise it
% would mean that this same code has to be copied for every quadrature type.
% This improve readability
function [nodes, weights] = calculate_quadrature(quad_function_handle, grid_level)
    nodes = cell(1, grid_level);
    weights = cell(1, grid_level);

    for i_order = 1:grid_level
        % Need to adjust the order, as a order of e.g. 3 would give 3 points,
        % but want then actually 5 points
        adjusted_order = 2*i_order - 1;
        [nodes{1, i_order}, weights{1, i_order}] = quad_function_handle(adjusted_order);
    end
end
