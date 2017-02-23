function pol_val = short_eval_pol(PCE, cubature)

    [u_pol_type, Ic_u_pol_type, I_u_pol_type] = unique(PCE.pol_type);

    % Get the unique quadratures
    temp_quads = cubature.quadratures.nodes(Ic_u_pol_type,:);

    temp_quads = unique(cell2mat(transpose(temp_quads)));

    %Determine values only for those quadratures
    temp_pol_val = evaluate_polynomial_u(u_pol_type, unique(PCE.basis(:,Ic_u_pol_type)), transpose(temp_quads));

    N_scen = size(cubature.scenarios,1);
    N_pols = size(PCE.basis,1);

    for i_u_pol_type = 1:numel(u_pol_type)
        cur_pol_order_index = I_u_pol_type == i_u_pol_type;
        N_cur_pols = nnz(cur_pol_order_index);

        [~, indexie1] = ismember(cubature.scenarios(:, cur_pol_order_index), temp_quads(:,i_u_pol_type));

       % To get the individual elements requires sub2ind, therefore need to
       % create all the elements that are needed
       pol_t_index = ones(N_cur_pols, N_pols, N_scen, 'int32');
       pol_index = repmat(int32(transpose(PCE.basis+1)), 1, 1, N_scen);    
       scen_index = repmat(permute(int32(indexie1), [2,3,1]), 1, N_pols, 1);

    %     indexie1 = repmat(permute(indexie1, [2, 3, 1]), 1, N_basis, 1);
        indexT = sub2ind_int(size(temp_pol_val), pol_t_index, pol_index, scen_index);

        pol_val = temp_pol_val(indexT);
    end
end
