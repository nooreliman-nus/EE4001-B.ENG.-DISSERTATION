function [ solution ] = IL_LP( interval, daily_price, duration, power, energy )
%This function aims to schedule non-interruptible loads
%   Inputs:
%       interval - how long a time interval is in minutes
%       daily_price - electricity prices(given in half hour intervals)
%       power - the power rating of the loads
%       energy - required energy of loads
%   Outputs:
%       F - objective function value
%       x - optimal schedule

interval = 10; %time step in minutes
daily_price = ((30/interval)/1000)*xlsread('ElectricityPrice.xlsx'); %$/kWh in half-hour resolution
daily_price = repelem(daily_price, 30/interval); %repeats each element 3 times to match 10-min resolution
N = length(daily_price); %no. of periods in 24 hours
power = 1;
energy = 10;

%Objective Function

f = daily_price * power; %cost

%Lower and upper bounds

LB = zeros(N,1);
UB = power * ones(N,1);

%   Equality Constraints
%       Aeq,beq are the matrix and vector components respectively

Aeq = ones(1,N);
beq = energy/interval;

%   Inequality Constraints
%       A,b are the matrix and vector components respectively

A = eye(N)
b = power * ones(N,1);

[solution] = linprog(f',A,b,Aeq,beq,LB,UB);
end

