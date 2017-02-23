function basis_norm = calculate_norm(basis, SettingsPCE)
%   CALCULATE_NORM calculates the norm of basis vectors
%
%   basis_norm = CALCULATE_NORM(basis, pol_order, pol_type) 
%   
%   INPUT: 
%          basis: the basis vectors for which to determine the norm [N x M
%          matrix, where N is the number of basis vectors and M the number
%          of dimensions]
%          pol_order: The maximum polynomial order for the basis [integer]
%          pol_type: The polynomial types of the different dimensions [M x
%          1 cell array]
%
%   OUTPUT: 
%          basis_norm: The norm of the basis vectors [N x 1 array]
    
    N_basis_vectors = size(basis, 1);
    pol_order = double(SettingsPCE.pol_order);
    
    N_pol_type = numel(SettingsPCE.pol_type);
    
    basis_norm = zeros(size(basis));
    
    % TODO: It is possible to do unique polynomials for everything that is
    % not arbitrary
    
    for i_pol_type = 1:N_pol_type
       % Determine the basis norm just for the current polynomial type
       % Norms are such that the sum of weights is 1
       switch SettingsPCE.pol_type{i_pol_type}
           case 'hermite'
               cur_basis_norm = factorial(0:pol_order);
           case 'legendre'
               cur_basis_norm = 1./(2.*(0:pol_order)+1);
           case 'laguerre'
               cur_basis_norm = ones(1, pol_order+1);
           case 'jacobi'
               pol_orders = 0:pol_order;
               numerator = 0.5.*2.*gamma(pol_orders+1).*gamma(pol_orders+1);
               denominator = (2.*pol_orders + 1).*factorial(pol_orders).*gamma(pol_orders +1);
               cur_basis_norm = numerator./denominator;
           case 'arbitrary'
               cur_basis_norm = transpose(SettingsPCE.pdf_settings.pdf_parameters{i_pol_type}.norm(1:pol_order+1));
          otherwise
               error('calculate_norm:unknownPolType','The requested polynomial type is not known');
       end % End switch over pol_types   
       
       % Determined for only one polynomial, need to expand to the actual number
       % of polynomials
       N_cur_pols = 1;
       cur_basis_norm = repmat(cur_basis_norm, N_cur_pols, 1);
       
       % Then find the positions of the data we actually need
       I_dim = repmat(1:N_cur_pols, N_basis_vectors, 1);
       I_pol = basis(:, i_pol_type) + 1;       
       I_cur_to_full = sub2ind(size(cur_basis_norm), I_dim, double(I_pol));

       basis_norm(:, i_pol_type) = cur_basis_norm(I_cur_to_full);

    end % End for loop over pol_types
        
    % Now need the product over the different polynomials
    basis_norm = prod(basis_norm, 2); 
end
