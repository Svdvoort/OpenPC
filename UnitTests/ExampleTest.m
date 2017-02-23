classdef ExampleTest < matlab.unittest.TestCase   
    % Unit tests for Gaussian quadrature
    methods (Test)
        function testExampleRun(testCase)
            settings = 'example_settings.json';
            PCE = OpenPC(settings);
            scenarios = [0,0,0,0; ...
                         1,1,1,1;
                         2,2,2,2];      
            PCE_result = evaluate_PCE(PCE, scenarios);
            % Need to round because of machine precision
            PCE_result = round(PCE_result, 10);
            blackbox_output = example_problem(scenarios, PCE.SettingsPCE.blackbox_arguments);
            blackbox_output = transpose(blackbox_output);

            testCase.verifyEqual(PCE_result, blackbox_output);
        end
    end
end 