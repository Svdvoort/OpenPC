function ExampleRun()
% This function provides an example run of constructing and evaluating a
% PCE using OpenPC

% Determine which settings to use. As long as the settings are in the
% settings folder, only need to provide the name, not the path of the file

settings = 'example_settings.json';

% Run OpenPC to construct the PCE
PCE = OpenPC(settings);

% The PCE is not yet automatically saved by OpenPC so if we want to use it
% later on we need to save it. Uncomment the following line to do that

%save('Example.mat', 'PCE');

% Now we can evaluate the constructed PCE, for the desired 'scenarios'
% This are the input values for the PCE function
% In this case we evaluate the PCE function for three scenarios:
% One where all 4 inputs are 0, where all 4 input are 1 and where all 4
% inputs are 2.
scenarios = [0,0,0,0; ...
             1,1,1,1;
             2,2,2,2];      
         
% Just provide the PCE and the scenarios to get the result
PCE_result = evaluate_PCE(PCE, scenarios);
% We will actually run into machine precision, so round for better display
if verLessThan('matlab', '8.4')
    PCE_result = round(PCE_result*1e10)/1e10;
else
    PCE_result = round(PCE_result, 1e10);
end

% Calculate the output of the real, blackbox function
blackbox_output = example_problem(scenarios, PCE.SettingsPCE.blackbox_arguments);

% Output and compare the results
for i = 1:3
    fprintf('The PCE result for scenario %d is: %d\n', i, PCE_result(i))
    fprintf('The real output for scenario %d is: %d\n', i, blackbox_output(i))
end

% Determine the difference between the real and PCE output
diff_PCE_real = abs(PCE_result - transpose(blackbox_output));

% If it's small enough, the ExampleRun successfully completed
if all(diff_PCE_real < 1e-8)
    disp('PCE results match the real function output!')
    disp('Example Run successfully completed!')
else
    % Output a warning, something is probably not right
    warning('PCE results did not match the real function output! Something might be wrong!')
end


end