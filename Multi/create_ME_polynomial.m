function SettingsPCE = create_ME_polynomial(SettingsPCE)      
%   CREATE_ME_POLYNOMIAL constructs the settings for multi-element
%   (arbitrary) PDFs
%
%   SettingsPCE = create_ME_polynomial(SettingsPCE)
%   
%   INPUT: 
%          SettingsPCE: The settings used to constructs the PCE [struct]
%   OUTPUT: 
%          SettingsPCE: Settings updated with multi-element settings
%          [struct]

    % Loop over the number of elements
    for i_element = 1:SettingsPCE(1).SettingsME.N_elements
        
        % The order for which we need to compute the multi-element settings
        % depend the polynomial order, as well as the grid order
        required_order = max(SettingsPCE(i_element).pol_order+1, 2*SettingsPCE(i_element).grid_level);

        % Within this element loop over all arbitrary pdfs
        for i_element_pdf = SettingsPCE(i_element).SettingsME.I_arbitrary

            % Beta, gamma and normalization factor are calculated using
            % stieltjes procedure
            [beta, gamma, norm] =...
                stieltjes(SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf},...
                                      required_order);

            % Update the settings with the calculated polynomials
            SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.beta = beta;
            SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.gamma = gamma;
            SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.norm = norm;
        end
    end
end