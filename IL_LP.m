function [ solution ] = IL_LP( t, pr, pw, e )
%This function aims to schedule interruptible loads without electrical
%machinery
%   Inputs:
%       t - how long a time interval is in minutes
%       pr - electricity prices(given in half hour intervals)
%       pw - the power rating of the loads
%       e - the energy required
%   Outputs:
%       F - objective function value
%       x - optimal schedule

N = length(pr); %no. of periods in 24 hours, no. of variables

%Objective Function - min f'*x

f = t * pr; %cost matrix

%   Equality Constraints
%       Aeq,beq are the matrix and vector components respectively

Aeq = ones(1,N); %array of ones corresponding to total load consumption = e
beq = e/t;

%   Inequality Constraints - A*x <= b
%       A,b are the matrix and vector components respectively
%       Include upper and lower bounds for each variable

A_bound_up = eye(N);
A_bound_low = -1* A_bound_up;

A = [A_bound_up;A_bound_low];

b_bound_up = pw * ones(N,1);
b_bound_low = zeros(N,1);

b = [b_bound_up;b_bound_low];

[solution,Final_Cost,exitflag] = linprog(f',A,b,Aeq,beq);

display(Final_Cost)
display(solution)
end

