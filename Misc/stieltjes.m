function [beta, gamma, norm] = stieltjes(pdf_settings, pol_order)
%   STIELTJES construct orthonomal polynomials using the Stieltjes
%   procedure
%
%   [beta, gamma, norm] = STIELTJES(pdf_function, pol_order, trunc_limit) returns
%   the alpha and beta of the constructed polynomial
%   
%   INPUT: 
%          pdf_function: Name of pdf function [string]
%          pol_order: The maximum polynomial order up to which to make the
%                     polynomial [Integer]
%          trunc_limit: The range of the current section (from lower to
%                     upper limit) [2x1 array]
%   OUTPUT: 
%          alpha: The alpha of the constructed polynomial [pol_order x 1
%                 array]
%          beta: The beta of the constructed polynomial [pol_order x1
%                array]
%          norm: The norm of the polynomial

    % TODO: These values are hardcoded
    accuracy = 1e-12;    
    N_points = 1000;
    
    pdf_function = str2func(pdf_settings.pdf_function);      

    % Determine beta, gamma and norm from integration using gauss_legendre
    % quadrature
    [quad_points, quad_weights] = gauss_legendre_quadrature(N_points, pdf_settings);
    quad_points = transpose(quad_points);
    
    % Calculate PDF value beforehand as well
    pdf_value = pdf_function(quad_points, pdf_settings);      
    
    % Adjust weights because going from -1 to 1 legendre_quadrature is
    % normalized in [-1, 1]
    quad_weights = quad_weights*(pdf_settings.b - pdf_settings.a);

    % Initialize variables
    norm = zeros(pol_order, 1);
    beta = zeros(pol_order, 1);
    gamma = zeros(pol_order, 1);
    function_value = ones(pol_order, N_points);

     for i_order = 1:pol_order        
        if i_order == 1
            % We already initialized ones, don't actually need to do
            % anything
        elseif i_order == 2
            function_value(i_order,:) = (quad_points + beta(i_order - 1)).*function_value(i_order - 1,:);
        else
            function_value(i_order,:) = (quad_points + beta(i_order - 1)).*function_value(i_order - 1,:)...
                                        - gamma(i_order - 1).*function_value(i_order - 2,:);
        end

        norm(i_order) = (function_value(i_order,:).^2.*pdf_value)*quad_weights;
        
        beta(i_order) = -(quad_points.*function_value(i_order,:).^2.*pdf_value)*quad_weights./norm(i_order);  
        
        if i_order == 1
            gamma(i_order) = norm(i_order);
        else
            gamma(i_order) = norm(i_order)/norm(i_order-1);        
        end
    end

    % All small element are set to 0
    beta(abs(beta) < accuracy) = 0;
    gamma(abs(gamma) < accuracy) = 0;
end