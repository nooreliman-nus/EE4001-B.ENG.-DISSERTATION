function [solution] = NL_LP(t, pr, l, pw)
%This function aims to schedule non-interruptible loads
%   Inputs:
%       t - how long a time interval is in minutes
%       pr - electricity prices(given in half hour intervals)
%       l - required time the loads have to run
%       pw - the power rating of the loads
%   Outputs:
%       F - objective function value
%       x - optimal schedule

N = length(pr); %no. of periods in 24 hours, no. of variables

%Objective Function - min f'*x

f = t * pr * pw; %cost matrix

%   Equality Constraints - Aeq*x = beq
%       Aeq,beq are the matrix and vector components respectively

Aeq = ones(1,N); %corresponds to total duration of ON load throughout l
beq = l;

%   Inequality Constraints - A*x <= b
%       A,b are the matrix and vector components respectively
%       Include upper and lower bounds for each variable

A_bound_up = eye(N);
A_bound_low = -1 * A_bound_up;

A = [A_bound_up;A_bound_low];

b_bound_up = ones(N,1);
b_bound_low = zeros(N,1);

b = [b_bound_up;b_bound_low];

%   Non-interruptibility Constraints
%       Create a band matrix where non-zero elements are in a diagonal band

first_col = [[l-1];[-l];zeros(N-l-1,1)];
first_row = [[l-1];-ones(l-1,1);zeros(N-l,1)];
band = toeplitz(first_col,first_row);

A = [A;band]; %Linear inequality constraint matrix

b = [b;zeros(N-l+1,1)]; %Linear inequality constraint vector

[solution,Final_Cost,exitflag] = linprog(f',A,b,Aeq,beq);

%Display final cost and power status
display(Final_Cost)
display(solution)

figure
%Plot of power status against time
subplot(2,1,1)
stairs(solution*pw,'black')
xlabel('Time')
ylabel('Power Status (W)')
%Plot of price against time
subplot(2,1,2)
plot(pr,'black')
xlabel('Time')
ylabel('Price ($/Wh)')

end

