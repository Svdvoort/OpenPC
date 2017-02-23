function [basis, norm] = create_basis_PC(PCE, SettingsPCE)
%   CREATE_BASIS_PC determines the basis vectors and norm of the basis vectors
%
%   [basis, norm] = CREATE_BASIS_PC(PCE, SettingsPCE) gives back the basis
%   and norm
%   
%   INPUT: 
%           PCE: structure with information about the PCE
%           [PCE structure]
%           SettingsPCE: Structure with settings for the construction of
%           the PCE [SettingsPCE structure]
%           
%
%   OUTPUT: 
%          basis: The different basis vectors to be included in the
%          expansion [N x M matrix, where N is number of basis vectors and
%          M is the number of dimension]
%          norm: Vector containing the norm of the basis vectors [N x 1
%          vector]

    N_pol_types = numel(PCE.pol_type);

    % The basis is basically just a sparse multi-index
    basis = create_sparse_multi_index(SettingsPCE.pol_order, uint8(N_pol_types));
    
    % If the user set trim, some basis vectors are automatically excluded such
    % that the coefficients of the remaining vectors can be accuretely determined
    if SettingsPCE.trim 
        N_dim = size(basis, 2);
        if N_dim == 1 
            trim_factor = 1;
        elseif SettingsPCE.grid_level < N_dim
            trim_factor = 1;
        else
            trim_factor = log(double(SettingsPCE.pol_order))./log(N_dim);
            trim_factor = trim_factor - 0.01;
        end

        % Here need to temporarily convert numeric type otherwise it is not
        % possible to raise to the power.
        I_to_keep = sum(single(basis).^trim_factor, 2).^(1/trim_factor) <= SettingsPCE.pol_order;
        basis = basis(I_to_keep,:);
        
    end

    % Also need to calculate the norm of the basis
    norm = calculate_norm(basis, SettingsPCE);
end
