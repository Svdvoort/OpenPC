function SettingsPCE = load_settings(settings_file)
%   LOAD_SETTINGS loads a json file with the settings for the PCE construction
%
%   SettingsPCE = load_settings(settings_file) returns the settings specified in
%   the settings file
%   
%   INPUT: 
%          settings_file: Name of settings file to use [String]
%   OUTPUT: 
%         SettingsPCE: Structure with the settings for PCE construction [Struct]

    SettingsPCE = loadjson(settings_file);
    
    % Settings load from json are strings, need to convert it to appropiate 
    % types    
    SettingsPCE.pol_order = uint8(SettingsPCE.pol_order);
    SettingsPCE.grid_level = uint8(SettingsPCE.grid_level);
    SettingsPCE.trim = logical(SettingsPCE.trim);
    SettingsPCE.remove_small_elements = logical(SettingsPCE.remove_small_elements);
    SettingsPCE.mem_saving = logical(SettingsPCE.mem_saving);
    
    SettingsPCE.pol_type = get_pol_type(SettingsPCE.input_distribution);
    SettingsPCE.quadrature_type = construct_quad_types(SettingsPCE.quadrature_type, SettingsPCE.pol_type);    
    
    % check to see if we're going to do multi-element (even if with only on element)    
    if any(ismember(SettingsPCE.input_distribution, 'arbitrary'))
        SettingsPCE.Do_ME = true;
        SettingsPCE.SettingsME = struct();
        SettingsPCE.SettingsME.I_arbitrary = find(ismember(SettingsPCE.input_distribution, 'arbitrary'));

        % See how many elements for each arbitrary pdf
        N_arbitrary = numel(SettingsPCE.SettingsME.I_arbitrary);
        I_element = cell(N_arbitrary, 1);
        for i_element_pdf = SettingsPCE.SettingsME.I_arbitrary       
           N_cur_element = numel(SettingsPCE.pdf_settings.pdf_parameters{i_element_pdf}.elements);
           I_element{i_element_pdf} = 1:N_cur_element;
        end
        
        % Get element from cartesian product of all indvidual elements       
        I_combined_elements = cartesian(I_element{:});
        SettingsPCE.SettingsME.N_elements = size(I_combined_elements, 1);
        SettingsPCE.SettingsME.N_arbitrary = N_arbitrary;

        full_pdf_parameters_arbitrary = {SettingsPCE.pdf_settings.pdf_parameters{SettingsPCE.SettingsME.I_arbitrary}};
        SettingsPCE = repmat(SettingsPCE, SettingsPCE.SettingsME.N_elements, 1);

        % Make the settings according to the 
        for i_element = 1:SettingsPCE(1).SettingsME.N_elements
          for i_element_pdf = 1:N_arbitrary
              i_cur_element = I_combined_elements(i_element, i_element_pdf);

               SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.a = ...
                   full_pdf_parameters_arbitrary{i_element_pdf}.elements{i_cur_element}.a;

               SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.b = ...
                   full_pdf_parameters_arbitrary{i_element_pdf}.elements{i_cur_element}.b;
               
               SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.elements = ...
                   [SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.a,...
                    SettingsPCE(i_element).pdf_settings.pdf_parameters{i_element_pdf}.b];

          end
        end
    else
        SettingsPCE.SettingsME = struct();
        SettingsPCE.SettingsME.N_elements = 1;
        SettingsPCE.Do_ME = false;
    end
end