function PCE_std = calculate_std(PCE)
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

    if PCE(1).SettingsPCE.Do_ME
        N_elements = numel(PCE);
        N_outputs = size(PCE(1).coeffs, 1);
        PCE_std = zeros(N_elements, N_outputs);
        for i_element = 1:N_elements
            PCE_std(i_element, :) = PCE(i_element).coeffs(:, 2:end)*PCE.norm(2:end);
        end
        
        % Need to apply a correction in case the arbitrary polynomial is
        % scaled
        for i_pdf_element = 1:PCE(1).SettingsPCE.SettingsME.N_arbitrary
            if isfield(PCE(1).SettingsPCE.pdf_settings.pdf_parameters{i_pdf_element}, 'std')
                PCE_std = PCE_std./PCE(1).SettingsPCE.pdf_settings.pdf_parameters{i_pdf_element}.std;
            end
        end

    else
        PCE_std = PCE.coeffs(:, 2:end)*PCE.norm(2:end);
    end
end