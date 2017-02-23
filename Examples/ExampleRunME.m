function ExampleRunME()
% This function provides an example run of constructing and evaluating a
% Multi-element PCE using OpenPC

% Determine which settings to use. As long as the settings are in the
% settings folder, only need to provide the name, not the path of the file

settings = 'example_settings_ME.json';

% Run OpenPC to construct the PCE
PCE = OpenPC(settings);

% The PCE is not yet automatically saved by OpenPC so if we want to use it
% later on we need to save it. Uncomment the following line to do that

%save('Example.mat', 'PCE');

% Now we can evaluate the constructed PCE, for the desired 'scenarios'
% This are the input values for the PCE function
% In this case we evaluate the PCE function for four scenarios;
% Three for the three different elements that are defined in the example
% problem and one for the case where it is outside of all elements, so it
% returns NaN
scenarios = [-5, -2, 0;
             0 ,1, 1;
             5, 2, 1;
             9, 5, 1];      
         
% Just provide the PCE and the scenarios to get the result
PCE_result = evaluate_PCE(PCE, scenarios);
% We will actually run into machine precision, so round for better display
if verLessThan('matlab', '8.4')
    PCE_result = round(PCE_result*1e10)/1e10;
else
    PCE_result = round(PCE_result, 1e10);
end

N_scenarios = size(scenarios, 1);
diff_PCE_real = ones(N_scenarios, 1);
% Output and compare the results
for i = 1:N_scenarios
    % Calculate the output of the real, blackbox function
    blackbox_output = example_problem_ME(scenarios(i, :), PCE(1).SettingsPCE.blackbox_arguments);
    
    % Determine the difference between the real and PCE output
    if i == 4
        if isnan(PCE_result(i))
            diff_PCE_real(i) = 0;
            fprintf('This element is not in the constructed PCE and should thus return NaN\n')
        end
    else
        diff_PCE_real(i) = abs(PCE_result(i) - transpose(blackbox_output));
    end

    fprintf('The PCE result for scenario %d is: %d\n', i, PCE_result(i))
    fprintf('The real output for scenario %d is: %d\n', i, blackbox_output)
    fprintf('\n')
end


% If it's small enough, the ExampleRun successfully completed
if all(diff_PCE_real < 1e-8)
    disp('PCE results match the real function output!')
    disp('Example Run successfully completed!')
else
    % Output a warning, something is probably not right
    warning('PCE results did not match the real function output! Something might be wrong!')
end


end