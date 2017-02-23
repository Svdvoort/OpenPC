function [result, temp_recur] = jacobi(input, pol_orders, pdf_settings)
%   HERMITE evaluates a probabilist' Hermite polynomial
%
%   result = HERMITE(input,pol_order) calculates the result of
%   evaluating a hermite polynomial with polynomial order pol_order at
%   the points input. 
%   
%   INPUT: 
%           input: Points at which to evaluate the Hermite polynomial
%           [N x M matrix, N is number of inputs, M is number of scenarios]
%           pol_order: The UNIQUE polynomial orders of the Hermite polynomial
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
        error('jacobi:noUniquePolynomials','This function only supports unique polynomial orders!')
    end

    N_input = size(input, 1);
    N_scenarios = size(input, 2);

    max_pol_order = max(pol_orders); 
    
    alpha = pdf_settings.alpha;
    beta = pdf_settings.beta;
    
    % Adjust input according to pdf_settings
    pdf_range = pdf_settings.b - pdf_settings.a;
    pdf_middle = (pdf_settings.b + pdf_settings.a);
    
    input = (input - pdf_middle)/pdf_range;

    
    % Need to reshape a bit to get it conform the shape of temp_recur
    input = permute(input, [1, 3, 2]);
    
    if max_pol_order == 0
        result = ones(N_input, 1, N_scenarios);
    else
        temp_recur = ones(N_input, max_pol_order+1, N_scenarios); 
        temp_recur(:, 2, :) = (alpha-beta)/2 + (1+ (alpha+beta)/2).*input;  
        
        for cur_pol = 3:max_pol_order+1
            n = cur_pol - 1;
            factor1 = (2*n + alpha + beta - 1).*((2*n + alpha + beta)*(2*n + alpha + beta - 2).*input + alpha^2 - beta^2);
            factor2 = 2.*(n + alpha - 1).*(n+beta-1).*(2*n + alpha + beta);
            factor3 = 2*n*(n + alpha + beta)*(2*n + alpha + beta -2);
            temp_recur(:, cur_pol, :) = (factor1.*temp_recur(:, cur_pol - 1, :)...
                - factor2.*temp_recur(:, cur_pol - 2, :))./factor3;
        end
        
        % Finally just return the values that are actually requested
        result = temp_recur(:, pol_orders + 1, :);        
    end
end
