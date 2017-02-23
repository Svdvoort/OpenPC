function output = evaluate_PCE(PCE, scenarios)
    %   EVALUATE_PCE evaluates a PCE for the given scenarios
    %
    %   output = EVALUATE_PCE(PCE, scenarios) calculates the result of
    %   evaluating the PCE given in the PCE structure for the supplied
    %   scenarios
    %   
    %   INPUT: 
    %           PCE: Structure containing all the information about the PCE
    %           Scenarios: The scenarios for which to sample the PCE [N x M
    %           matrix, where N is the number of scenarios and M the number
    %           of inputs in the PCE]
    %
    %   OUTPUT: 
    %           Output: The result of evaluating the PCE. Return Nan when 
    %           doing multi_element and the range is not in one of the 
    %           elements [N x K matrix, K is number of responses in PCE]

    % First we check whether we have multi-element PCE, because then we
    % have to select which element we're actually going to use
    N_elements = numel(PCE);
    if N_elements > 1
        output = NaN(size(scenarios, 1), size(PCE(1).coeffs, 1));
        I_arbitrary = PCE(1).SettingsPCE.SettingsME.I_arbitrary;
        arbitrary_inputs = scenarios(:, I_arbitrary);
        for i_element = 1:N_elements
            element_limits = ...
                cell2mat(cellfun(@(c) [c.a; c.b], {PCE(i_element).SettingsPCE.pdf_settings.pdf_parameters{I_arbitrary}},...
                'Uniform', 0));
            
            in_range = bsxfun(@ge, arbitrary_inputs, element_limits(1,:)) &...
                           bsxfun(@lt, arbitrary_inputs, element_limits(2,:));
            in_range = logical(prod(in_range, 2));
            
            in_range_scenarios = scenarios(in_range, :);
            
            if numel(in_range_scenarios) > 0
                pol_val = evaluate_polynomial(PCE(i_element).pol_type, PCE(i_element).basis, transpose(in_range_scenarios), PCE(i_element).SettingsPCE);
                output(in_range, :) = PCE(i_element).coeffs*transpose(pol_val);
            end
        end
    else
        % Evaluate the polynomials themselves
        pol_val = evaluate_polynomial(PCE.pol_type, PCE.basis, transpose(scenarios), PCE.SettingsPCE);
        % Finally the output is just the value of the polynomial times the PCE
        % coefficients
        output = PCE.coeffs*transpose(pol_val);
    end
        

    
end
