function output = example_problem_ME(input, blackbox_arguments)
%  EXAMPLE_PROBLEM_ME is an example interface of the blackbox
%
%   output = EXAMPLE_PROBLEM_ME(input, optional_arguments)
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
    switch blackbox_arguments.problem_number
        case 1 
            % A very simple polynomial. The second index is the different inputs.
            if all(input(:, 1) <= -2)
                output = ones(size(input, 1), 1);
            elseif all(input(:,1) > -2) && all(input(:,1)< 2)
                output = 2*input(:, 1);
            else
                output = input(:, 1);
            end
        case 2
            if all(input(:,1) <= -3)
                output = 2.*input(:,1);
            elseif all(input(:,1) >= -3) &&all(input(:, 1) <= 3)
                output = input(:,1).^2;
            else
                output = 2.*input(:,1).^3;
            end
        case 3
            % Problem with multiple multi-element inputs
            if all(input(:,1) <= -3) && all(input(:, 2) <= -1)
                output = 2.*input(:,1) + input(:,2).^2;
            elseif all(input(:,1) >= -3) && all(input(:, 1) <= 3)
                output = input(:,1).^2;
            else
                output = 2.*input(:,3).^3;
            end
    end
end