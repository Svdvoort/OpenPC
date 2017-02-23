function SettingsPCE = create_arbitrary_polynomial(SettingsPCE)     

    for i_element = 1:SettingsPCE(1).SettingsME.N_elements
        required_order = max(SettingsPCE(i_element).pol_order, 2*SettingsPCE(i_element).grid_level);
        
        trunc_limit = (SettingsPCE(i_element).SettingsME.elements - SettingsPCE.means)*SettingsPCE.stds;

        [beta, gamma, norm] = stieltjes(SettingsPCE(i_element).SettingsME.pdf_function, required_order, ...
                                    trunc_limit,...
                                    SettingsPCE(i_element).SettingsME.pdf_settings);

        SettingsPCE(i_element).SettingsME.beta = beta;
        SettingsPCE(i_element).SettingsME.gamma = gamma;
        SettingsPCE(i_element).SettingsME.norm = norm;
    end
            

end