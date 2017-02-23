function multi_index = create_sparse_multi_index(max_value, N_dimensions)
%   CREATE_SPARESE_MULTI_INDEX createes a sparse multi index
%
%   multi_index = CREATE_SPARE_MULTI_INDEX(max_value, N_dimensions) creates
%   a matrix such that the sum of each row <= max_value, where a row
%   has a number of elements equal to N_dimensions
%   
%   INPUT: 
%           max_value: Maximum value for sum of each row [Integer]
%           N_dimension: Number of elements for each row [Integer]
%
%   OUTPUT: 
%           multi_index: The constructed multi-index [NxM matrix where N =
%           number of combinations and M = N_dimensions] 

    % Initialize cell for construction of basis, this is done progressively
    temp_multi_index = cell(max_value + 1, N_dimensions);

    % When only 1 dimension it's trivial, just up to polynomial order
    temp_multi_index(:,1) = num2cell(0:max_value);

    for i_dim = 2:N_dimensions % Loop over all of the dimensions
        for i_value = 0:max_value % For each dimension, loop up to the maximum value
            for i_intermediate_value = 0:i_value
                % Get the previously constructed 
                temp_result = temp_multi_index{i_value - i_intermediate_value + 1, i_dim-1}; 
                temp_result(:, end + 1) = i_intermediate_value; %#ok add the current value
                % Update the  matrix with new values
                temp_multi_index{i_value + 1, i_dim} = [temp_multi_index{i_value + 1, i_dim}; temp_result]; 
            end        
        end
    end

    % Need the final result as a matrix
    multi_index = cell2mat(temp_multi_index(:, N_dimensions));
end