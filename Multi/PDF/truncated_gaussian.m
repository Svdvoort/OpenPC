function pdf_values = truncated_gaussian(x, pdf_settings)

    std = pdf_settings.std;
    mean = pdf_settings.mean;
    
    lower_lim = pdf_settings.a;
    upper_lim = pdf_settings.b;
    
    rescaled_pdf = 1/std.*normpdf((x - mean)./std);

    % Need to normalize to get total area of 1
    normalization_factor = normcdf((upper_lim - mean)/std)...
                           - normcdf((lower_lim - mean)/std);

    pdf_values = rescaled_pdf./normalization_factor;
    
    % All values outside the current range are set to 0
    outside_range = x <= lower_lim | x >= upper_lim;
    pdf_values(outside_range) = 0;
end