function sparse_grid = make_sparse_grid(cubature)

    sparse_grid.nodes{1} = cubature.nodes{1};
    sparse_grid.weights{1} = cubature.weights{1};
    
    for i = 2:length(cubature.nodes)
       sparse_grid.nodes{i} = union(cubature.nodes{i-1}, cubature.nodes{i}, 'rows');
       
       weightTempA = ismember(sparse_grid.nodes{i}, cubature.nodes{i-1}, 'rows');
       weightTempB = ismember(sparse_grid.nodes{i}, cubature.nodes{i}, 'rows' );
       
       sparse_grid.weights{i} = zeros(size(sparse_grid.nodes{i},1),1);
       
       sparse_grid.weights{i}(weightTempA) = sparse_grid.weights{i}(weightTempA) - cubature.weights{i-1};
       sparse_grid.weights{i}(weightTempB) = sparse_grid.weights{i}(weightTempB) + cubature.weights{i};

    end

end
