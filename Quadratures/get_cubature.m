function cubature = get_cubature(quadrature, grid_level, grid_type)
%   GET_CUBATURE determines the nodes and weights for a cubature
%
%   cubature = GET_CUBATURE(quadrature, grid_level, grid_type) gives back the
%   cubature structure with the nodes, weights and scenarios
%   
%   INPUT: 
%           quadrature: Contains the nodes and weights of the quadratures
%           [quadrature structure]
%           grid_level: The maximum grid level [Integer]
%           grid_type: The type of grid to be used. Supported are:
%           'full','sparse' and 'sparse-adaptive'
%
%   OUTPUT: 
%          cubature: Contains the nodes, weights and scenarios of the cubature
%          [Cubature structure]
    
    cubature = struct();
    cubature.grid_type = grid_type; 
    cubature.quad_rules = quadrature.quad_rules;
    
    % First make the multi_index
    % Basically just the combinations of the different polynomial orders
    N_quad_types = uint8(size(quadrature.nodes, 1));
    % Now different grid types are possible, this influences the multi index
    % that we keep
    switch grid_type
        case 'full'
            % In case of the full grid, just the need the highest multi_index,
            % as this contains everything
            multi_index = cartesian(repmat(1:grid_level, uint8(N_quad_types), 1));   
            multi_index_norm = sum(multi_index, 2, 'native');
            to_keep = multi_index_norm == max(multi_index_norm);
            multi_index = multi_index(to_keep, :);
        case 'sparse'
            % With sparse we can create sparse grids, don't need to compute
            % full grid first. We need to add 1 to the constructed multi
            % index because we start at grid 1 instead of 0
            multi_index = create_sparse_multi_index(grid_level-1,N_quad_types)+uint8(1);
        case 'sparse-adaptive'
            % For future use
            multi_index = create_sparse_multi_index(grid_level-1,N_quad_types)+uint8(1);
        otherwise
            error('get_cubature:UnknownGridType','The requested grid type is not known!');
    end
     
    cubature.multi_index = multi_index;

    
    % Now determine which of the quadratures need to be combined based on
    % the multi index we just obtained
    % The quadratures that are needed that easily be obtained, just need to
    % calculate the correct index, need dimension index for that (which
    % just says that input 1 is at column 1, input 2 at column 2 etc)
    I_dim = repmat(1:N_quad_types, size(multi_index, 1), 1);    
    needed_quad = transpose(sub2ind(size(quadrature.nodes), I_dim, multi_index));
    
    cubature.nodes = cell(1, size(needed_quad,2));
    cubature.weights = cell(1, size(needed_quad,2));
    
    % Finally need to combine them, for the weight only need product
    for i_node = 1:size(needed_quad, 2)        
         cubature.nodes{1, i_node} = cartesian(...
             quadrature.sparse_nodes{needed_quad(:, i_node)});   
         cubature.weights{1, i_node} = prod(cartesian(...
             quadrature.sparse_weights{needed_quad(:, i_node)}), 2);
    end  
end
