function cubature = get_error_scenarios(SettingsPCE)
%   GET_ERROR_SCENARIOS determines the error scenarios used in PCE construction
%
%   cubature = GET_ERROR_SCENARIOS(SettingsPCE)
%   
%   INPUT: 
%           SettinsgPCE: Contains the different settings to construct the
%           error scenarios [SettingsPCE struct]
%
%   OUTPUT: 
%          cubature: Contains the nodes, weights and scenarios of the cubature
%          [Cubature structure]

    % First get the individual quadratures
    quadratures = get_quadratures(SettingsPCE.grid_level, SettingsPCE.quadrature_type, SettingsPCE.pdf_settings);
    
    % Now from the quadrature need to make cubature
    cubature = get_cubature(quadratures, SettingsPCE.grid_level, SettingsPCE.grid_type);
    cubature.quadratures = quadratures;
            
    % The cubature.nodes gives the scenarios for the different grids,
    % easier to have them in one big matrix, to calculate them all at once
    % instead of per grid.
    [cubature.scenarios, ~, scenario_index] = unique(cell2mat(transpose(cubature.nodes)), 'rows');
    
    % Instantely get the total weights for all scenarios, easier later on
    % instead of having to loop over grids.
    cubature.total_weights = accumarray(scenario_index, cell2mat(transpose(cubature.weights)));
    
    % Only need scenarios with non-zero weights for the calculation
    LI_non_zero_weights = cubature.total_weights ~= 0;
    
    cubature.scenarios = cubature.scenarios(LI_non_zero_weights, :);
    cubature.total_weights = cubature.total_weights(LI_non_zero_weights);
    
    % The scenarios will be parsed all at once to the black box, need to
    % keep track of where the different scenarios are in the grids.    
    node_size = cellfun('size', cubature.nodes, 1);
    cubature.full_index = mat2cell(scenario_index, node_size, 1);
end
