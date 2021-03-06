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


        % Calculate the coefficients
        if SettingsPCE.mem_saving
            % If memory saving is enabled, iterate over all scenarios
            % Smaller matrix sizes, but much slower
            coefficients = zeros(1, size(PCE.basis, 1));
        
            scenarios = transpose(cubature.scenarios);
           for i_scenario = 1:size(scenarios, 2)
             pol_val = evaluate_polynomial(PCE.pol_type, PCE.basis, scenarios(:, i_scenario), SettingsPCE);
             coefficients = coefficients + output(i_scenario, :)*pol_val*cubature.total_weights(i_scenario);
           end
        else
            % If not memory saving, fast calculation of coefficients
            % Using complete marix at once
            pol_val = evaluate_polynomial(PCE.pol_type, PCE.basis, transpose(cubature.scenarios), SettingsPCE);
            coefficients = transpose(output)*bsxfun(@times, pol_val, cubature.total_weights);
        end       
        
       
        coefficients = bsxfun(@rdivide, coefficients, transpose(PCE.norm));
        
        if SettingsPCE.remove_small_elements
            coefficients(abs(coefficients) < SettingsPCE.small_element_threshold) = 0;
        end
end
