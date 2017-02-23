function polynomial_value = evaluate_polynomial_new(pol_types, pol_orders, scenarios, SettingsPCE)
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
%          polynomials [K x M]
    
    N_pol_orders = size(pol_orders, 1);
    N_scen = size(scenarios, 2);
    
    % First find all cells with the same polynomial type, otherewise might
    % be doing unnecessary calculations
%     [u_pol_types, ~, I_u_pol_type] = unique(pol_types);      
%     N_u_pol_types = numel(u_pol_types);

    N_pol_types = numel(pol_types);    

    % Some pre-allocation
    polynomial_value = zeros(N_pol_orders, N_scen, N_pol_types);
    
    % Now iterate over all the unique polynomial types
    for i_pol_type = 1:N_pol_types
       % We need to gather the information for the current pol_type
%        LI_cur_pol_order = I_u_pol_type == i_u_pol_type;
%        N_cur_pols = nnz(LI_cur_pol_order);
       
%        cur_pol_orders = pol_orders(:, LI_cur_pol_order);
       U_cur_pol_orders = unique(pol_orders(:, i_pol_type));

       % Get only the unique polynomial orders, reduces calculations that
       % need to be done 
%        u_cur_pol_orders = double(unique(cur_pol_orders));
       
       % Now we evaluate the actual polynomial
       switch pol_types{i_pol_type}
           case 'hermite'
               polynomial_output = hermite(scenarios(i_pol_type,:), U_cur_pol_orders);
           case 'legendre'
               polynomial_output = legendre(scenarios(i_pol_type,:), U_cur_pol_orders);
           case 'laguerre'
               polynomial_output = laguerre(scenarios(i_pol_type,:), U_cur_pol_orders, SettingsPCE.pdf_settings.pdf_parameters{i_pol_type}.alpha);
           case 'jacobi'
               polynomial_output = jacobi(scenarios(i_pol_type,:), U_cur_pol_orders,  SettingsPCE.pdf_settings.pdf_parameters{i_pol_type}.alpha,  SettingsPCE.pdf_settings.pdf_parameters{i_pol_type}.beta);
           case 'arbitrary'
               polynomial_output = arbitrary_polynomial(scenarios(i_pol_type,:), U_cur_pol_orders, SettingsPCE.SettingsME.beta, SettingsPCE.SettingsME.gamma);
               
           otherwise
               error('evaluate_polynomial:unknownPolType','The requested polynomial type is not known');
       end     

       pol_output_size = size(polynomial_output);
       
       % To get the individual elements requires sub2ind, therefore need to
       % create all the elements that are needed
       I_pol_type = repmat(transpose(1:int32(N_cur_pols)), 1, N_pol_orders, N_scen);
       I_pol_order= repmat(transpose(int32(U_cur_pol_orders + 1)), 1, 1, N_scen);
       I_scenarios = repmat(permute(1:int32(N_scen),[3, 1, 2]), N_cur_pols, N_pol_orders, 1);
              
       I_cur_to_full = sub2ind_int(pol_output_size, I_pol_type, I_pol_order, I_scenarios);
       
       % Only need to store product for each polynomial type
       polynomial_value(:, :, i_u_pol_type) = prod(polynomial_output(I_cur_to_full), 1);
    end  % End of iteration over polynomial types
    % Take product over all polynomial types combined
    polynomial_value = transpose(prod(polynomial_value, 3)); 
end
