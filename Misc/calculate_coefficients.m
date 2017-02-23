function coefficients = calculate_coefficients(cubature, output, PCE, SettingsPCE)
%   CALCULATE_COEFFIENTS calculations the expansion coefficients of a PCE
%
%   coefficients = calculate_coefficients(cubature, output, PCE)
%   
%   INPUT: 
%           cubature: the cubatures to be used [cubature structure]
%           output: the blackbox model evaluated at the cubature nodes [N x
%           M matrix, where N is the number of nodes and M is the number of outputs]
%           PCE: contains information about the PCE [PCE structure]
%
%   OUTPUT: 
%          coefficients: The calculated coefficients of the expansion [M x
%          K matrix, where K is the number of basis vectors]

        % First calculate the polynomial value for all scenarios
        pol_val = evaluate_polynomial(PCE.pol_type, PCE.basis, transpose(cubature.scenarios), SettingsPCE);

        % Calculate the coefficients
        coefficients = transpose(output)*bsxfun(@times, pol_val, cubature.total_weights);
        coefficients = bsxfun(@rdivide, coefficients, transpose(PCE.norm));
        
        if SettingsPCE.remove_small_elements
            coefficients(abs(coefficients) < SettingsPCE.small_element_threshold) = 0;
        end
end
