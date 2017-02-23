function PCE_mean = calculate_mean(PCE)
    %   CALCULATE_MEAN extracts the mean from the PCE coefficients
    %
    %   PCE_mean = CALCULATE_MEAN(PCE) calculates the mean for the supplied
    %   PCE
    %   
    %   INPUT: 
    %           PCE: Structure containing all the information about the
    %           PCE [Structure]
    %
    %   OUTPUT: 
    %           PCE_mean: The mean of the output of the PCE 
    %           [N x K matrix, N is number of elements, 
    %            K is number of responses in PCE]

    N_elements = numel(PCE);
    if N_elements > 1
        N_outputs = size(PCE(1).coeffs, 1);
        PCE_mean = zeros(N_elements, N_outputs);
        for i_element = 1:N_elements
            PCE_mean(i_element, :) = PCE(i_element).coeffs(:, 1);
        end
    else
        PCE_mean = PCE.coeffs(:, 1);
    end
end