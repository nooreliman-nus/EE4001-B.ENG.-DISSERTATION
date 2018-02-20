function [solution] = TCL_LP(t, pr, pw, c, m, temp_up, temp_0, temp_req, temp_en, d)
%This function aims to schedule thermostatically controlled loads
%   Inputs:
%       t - how long a time interval is in minutes
%       pr - electricity prices(given in half hour intervals)
%       pw - the power rating of the loads
%       c- the specific heat of water
%       m - the mass of water in storage
%       temp_up - the upper limit of water temperature in storage
%       temp_0 - the initial water temperature in storage
%       temp_req - the required water temperature
%       temp_en - the environmental temperature at every interval
%       d - demand of hot water drawn at every interval
%   Outputs:
%       fval - objective function value
%       x - optimal schedule

N = length(pr); %no. of periods in 24 hours, no. of variables

%Objective Function - min f'*x

f = t * pr; %cost matrix

%   Inequality Constraints - A*x <= b
%       A,b are the matrix and vector components respectively
%       Include upper and lower bounds for each variable

%       x_tmp - to get the first i-s accumulated energy
    x_tmp = t * ones(N);
    x_tmp = tril(x_tmp);
%       x_tmp2 - opposite sign of x_tmp
    x_tmp2 = -1 * x_tmp;
    A = [x_tmp; x_tmp2];
%       ub - to obtain the power value
    ub = eye(N);
%       lb - opposite sign of ub
    lb = -1 * ub;
%       A_ub - coefficient in inequality A * x <= b
    A = [A;ub;lb];

%       C - heat consumption at each time step
    C = zeros(N); %initialise matrix)
    for i = 1:N
        C(i) = d(i) * c * (temp_req - temp_en(i);
    end
    
    b = zeros(2 * N, 1);
    for i = 1:N
%b(i, 0) - lower bound of energy accumulated at each time step
        b(i, 0) = sum(C(0:i + 1)) + m * c * (temp_up - temp_0);
%b(i + N, 0) - upper bound of energy accumulated at each time step
        b(i + N, 0) = -1 * sum(C(0:i + 1));
    end
    
%b_tmp - upper bound and then lower bound for power value at each time step
    b_tmp = zeros(2 * N, 1);
    b_tmp(0:N) = pw;
%b - value in inequality A * x <= b
    b = [b;b_tmp];
        
[x,fval,exitflag] = linprog(f',A,b,Aeq,beq);

end