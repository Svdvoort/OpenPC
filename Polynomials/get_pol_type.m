function pol_type = get_pol_type(input_distribution)
%   GET_POL_TYPE finds the polynomial types associated with the input
%   distribution
%
%   pol_type = GET_POL_TYPE(input_distribution) will return the polynomial
%   types for the input distributions given in input_distribution according
%   to the Wiener-Askey scheme
%   
%   INPUT: 
%           input_distribution: The input distributions for the different
%           inputs [1 x N cell array]
%
%   OUTPUT: 
%           pol_type: The associated polynomial types [1 X N cell array]

    N_inputs = numel(input_distribution);
    pol_type = cell(1, N_inputs);

    for i_input = 1:N_inputs
        switch input_distribution{i_input}
            case 'gaussian'
                pol_type{i_input} = 'hermite';
            case 'uniform'
                pol_type{i_input} = 'legendre';
            case 'gamma'
                pol_type{i_input} = 'laguerre'; 
            case 'beta'
                pol_type{i_input}  = 'jacobi';
            case 'arbitrary'
                pol_type{i_input} = 'arbitrary';
            otherwise
                error('get_pol_type:unknownInputDistribution','The input distribution does not have a polynomial type associated with it');

        end
    end
end