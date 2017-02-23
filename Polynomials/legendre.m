function [result, temp_recur] = legendre(input, pol_orders, pdf_settings)
%   LEGENDRE evaluates a Legendre polynomial
%
%   result = LEGENDRE(input,pol_order) calculates the result of
%   evaluating a Legendre polynomial with polynomial order pol_order at
%   the points input. 
%   
%   INPUT: 
%           input: Points at which to evaluate the Legendre polynomial
%           [N x M matrix, N is number of inputs, M is number of scenarios]
%           pol_order: The UNIQUE polynomial orders of the Legendre polynomial
%           [K x 1 vector, where K is the different polynomial orders]
%
%   OUTPUT: 
%           result: Result of evaluating the input for polynomial order
%           pol_order [length(input) X 1 vector]
%           temp_recur: Gives back not only the result for pol_order but
%           also for all lower orders [length(input) X pol_order + 1 matrix]

    % Let's do a quick sanity check to see if the given polynomial orders
    % are unique

    if numel(pol_orders) ~= numel(unique(pol_orders))
        error('legendre:noUniquePolynomials','This function only supports unique polynomial orders!')
    end

    N_input = size(input, 1);
    N_scenarios = size(input, 2);

    max_pol_order = max(pol_orders); 
    
    % Adjust input according to pdf_settings
    pdf_range = (pdf_settings.b - pdf_settings.a)/2;
    pdf_middle = (pdf_settings.b + pdf_settings.a)/2;
    
    input = (input - pdf_middle)/pdf_range;

    
    % Need to reshape a bit to get it conform the shape of temp_recur
    input = permute(input, [1, 3, 2]);
    
    if max_pol_order == 0
        result = ones(N_input, 1, N_scenarios);
    else
        temp_recur = ones(N_input, max_pol_order+1, N_scenarios); 
        temp_recur(:, 2, :) = input;  
        
        for cur_pol = 3:max_pol_order+1
            temp_recur(:, cur_pol, :) = ((2 .* cur_pol - 3) .* input .* temp_recur(:, cur_pol-1, :) - (cur_pol - 2) .* temp_recur(:, cur_pol-2, :)) ./ (cur_pol-1);
        end
        
        % Finally just return the values that are actually requested
        result = temp_recur(:, pol_orders + 1, :);        
    end
end