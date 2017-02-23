function output = example_problem(input, blackbox_arguments)
%  EXAMPLE_PROBLEM is an example interface of the blackbox
%
%   output = EXAMPLE_PROBLEM(input, optional_arguments)
%   
%   INPUT: 
%          input: The points for which the blackbox should be evaluated [N
%          x M matrix, where N is the number of scenarios and M the number
%          of inputs]
%          blackbox_arguments: Additional arguments which can be passed to
%          the function [struct]
%
%   OUTPUT: 
%          output: The output of the blackbox model [N x K matrix, where N is
%          the number of scenarios and K is the number of outputs]

    % Multiple problems are defined here
    % Very simple polynomials. The second index is the different inputs.

    switch blackbox_arguments.problem_number
        case 1 
            output = input(:, 1);
        case 2
            output = 2*input(:, 1).^2 + input(:,2);
        case 3
            output = 3*input(:,1) + 2*input(:,2).^3 + 3*input(:,3) + 21;
        case 4
            output = input(:,1).^2 + 3*input(:,2) + 2*input(:,3) +input(:,4).^2;

    end
end