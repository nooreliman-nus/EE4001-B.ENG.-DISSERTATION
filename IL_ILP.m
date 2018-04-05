function [ solution ] = IL_ILP( t, pr, pw, e )
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

A = [A_bound_up;A_bound_low]; %Linear inequality constraint matrix

b_bound_up = pw * ones(N,1);
b_bound_low = zeros(N,1);

b = [b_bound_up;b_bound_low]; %Linear inequality constraint vector

intcon = []; %There are no integer variables

[solution,Final_Cost,~] = intlinprog(f',intcon,A,b,Aeq,beq);

%Display final cost and power status
display(Final_Cost)
display(solution)

figure
%Plot of power status against time
subplot(2,1,1)
stairs(solution,'black')
xlabel('Time')
ylabel('Power Status (W)')
%Plot of price against time
subplot(2,1,2)
plot(pr,'black')
xlabel('Time')
ylabel('Price ($/Wh)')
end

