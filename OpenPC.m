function PCE = OpenPC(settings_file)
%   OpenPC creates a polynomial chaos expansion of a problem
%
%   PCE = OpenPC(settings_file) construct the PCE based on the given
%   settings
%   
%   INPUT: 
%          settings_file: Name of settings file to use [String]
%   OUTPUT: 
%          PCE: The constructed polynomial chaos expansion [PCE struct]

if verLessThan('matlab', '8.1')
    warning('OpenPC:Version', 'OpenPC has not been tested on Matlab versions before 2013a, correct functionality can not be guaranteed')
end

    % Initialize the OpenPC enviroment
    initOpenPC()

    SettingsPCE = load_settings(settings_file);
    
    % Initialize the PCE structure
    PCE = struct();
    PCE.pol_type = SettingsPCE(1).pol_type;
    PCE.quadrature_type = SettingsPCE(1).quadrature_type;    
    
    if SettingsPCE(1).Do_ME
        % Multi-element PCE, need to compute polynomials
        SettingsPCE = create_ME_polynomial(SettingsPCE);
        PCE = repmat(PCE, SettingsPCE(1).SettingsME.N_elements, 1);
    end
    
    for i_element = 1:SettingsPCE(1).SettingsME.N_elements
        % First add the basis to the structure
        [PCE(i_element).basis, PCE(i_element).norm] = create_basis_PC(PCE(i_element), SettingsPCE(i_element));          
        
        % First get all the information about the cubature
        cubature = get_error_scenarios(SettingsPCE(i_element));

        % Calculate the response of the blackbox
        output = blackbox(SettingsPCE(i_element).blackbox_function, cubature.scenarios, SettingsPCE(i_element).blackbox_arguments);    

        % Determine the coefficients    
        PCE(i_element).coeffs = calculate_coefficients(cubature, output, PCE(i_element), SettingsPCE(i_element));
        PCE(i_element).SettingsPCE = SettingsPCE(i_element);
    end 
end
