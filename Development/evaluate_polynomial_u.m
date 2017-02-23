function polynomial_value = evaluate_polynomial_u(pol_types, pol_orders, scenarios, varargin)
    %  EVALUATE_POLYNOMIAL evaluates a polynomial
    %
    %   polynomial_value = EVALUATE_POLYNOMIAL(pol_types, pol_orders, input)
    %   evaluates the polynomials given in pol_types for the polynomial
    %   orders given in pol_orders for the scenarios in input
    %   
    %   INPUT: 
    %           pol_types: Indicate the polynomial types for the different
    %           inputs [N x 1 cell array, N is number of inputs]
    %           pol_orders: Indicates the different polynomial orders to
    %           use [M x N matrix, where M is the total number of
    %           polynomials to evaluate]
    %           scenarios: Indicates the scenarios for which to evalute the
    %           polynomials [N x K matrix, where K is the number of
    %           scenarios]
    %
    %   OUTPUT: 
    %          polynomial_value: Gives the result of evaluating the
    %          polynomials [N x M x K]
    
    %   This software is released under the <a href="matlab: 
    %   web('https://www.gnu.org/licenses/gpl-3.0.en.html')">GPLv3 license</a>.
    %   Copyright (C) 2015  Sebastian van der Voort
    
    % First find all cells with the same polynomial type, otherewise might
    % be doing unnecessary calculations
    
    N_pols_t = numel(pol_types);
    N_pols = size(pol_orders, 1);
    N_scen = size(scenarios, 2);

    % Some pre-allocation
    polynomial_value = zeros(N_pols_t, N_pols, N_scen);
    
    % Now iterate over all the unique polynomial types
    for i_pol_type = 1:N_pols_t
       % We need to gather the information for the current pol_type
       cur_pol_order_index = i_pol_type;
       N_cur_pols = nnz(cur_pol_order_index);

       cur_pol_orders = pol_orders(:, cur_pol_order_index);
       % Get only the unique polynomial orders, reduces calculations that
       % need to be done 
       u_cur_pol_orders = double(cur_pol_orders);
       
       % Now we evaluate the actual polynomial
       switch pol_types{i_pol_type}
           case 'hermite'
               polynomial_output = hermite(scenarios, u_cur_pol_orders);
           case 'legendre'
               polynomial_output = legendre(scenarios, u_cur_pol_orders);
           case 'laguarre'
               polynomial_output = laguerre(scenarios, u_cur_pol_orders);
           otherwise
               error('evaluate_polynomial:unknownPolType','The requested polynomial type is not known');
       end     

       pol_output_size = size(polynomial_output);
       
       % To get the individual elements requires sub2ind, therefore need to
       % create all the elements that are needed
       pol_t_index = repmat(transpose(1:int32(N_cur_pols)), 1, N_pols, N_scen);
       pol_index = repmat(transpose(int32(cur_pol_orders + 1)), 1, 1, N_scen);
       scen_index = repmat(permute(1:int32(N_scen),[3, 1, 2]), N_cur_pols, N_pols, 1);
              
       indexT = sub2ind_int(pol_output_size, pol_t_index, pol_index, scen_index);
       
       polynomial_value(cur_pol_order_index, :, :) = polynomial_output(indexT);
    end  % End of iteration over polynomial types
end
