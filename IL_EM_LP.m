function [ solution ] = IL_EM_LP( t, pr, pw, e, t_off, theta )
%This function aims to schedule interruptible loads without electrical
%machinery
%   Inputs:
%       t - how long a time interval is in minutes
%       pr - electricity prices(given in half hour intervals)
%       pw - the power rating of the loads
%       e - the energy required
%       t_off - the minimum off time
%       theta - a small constant to connect z variables
%   Outputs:
%       F - objective function value
%       x - optimal schedule
%       y - optimal schedule of ancilliary variables
%       z - optimal schedule of constraint variables

N = length(pr); %no. of periods in 24 hours, no. of variables

%Objective Function - min f'*x

f_x = t * pr; %corresponds to power status x
f_y = zeros(1, N); %corresponds to ancilliary binary variables y
f_z = zeros(1, 3*N); %corresponds to variable z connecting x and y
f = [f_x, f_y ,f_z]; %concatenates f_x, f_y and f_z to form the objective function

%   Equality Constraints
%       Aeq,beq are the matrix and vector components respectively

Aeq = [ones(1, N),zeros(1, 4 * N)];
%y_tmp - sub-vector of Aeq to obtain y value at each time step
    y_tmp = eye(N);
%z_tmp1 - sub-vector of Aeq to obtain -(z2 + z3) at each time step
    z_tmp1 = zeros(N, 3 * N);
    for i = 1:N
        z_tmp1(i, ((i-1) * 3 + 2): ((i) * 3)) = -1;
    end
%tmp1 - coefficient to satisfy equality constraint: y - z2 - z3 = 0
    tmp1 = [zeros(N), y_tmp, z_tmp1];
%x_tmp - sub-vector of Aeq to obtain x value at each time step
    x_tmp = eye(N);
%z_tmp2 - sub-vector of Aeq to obtain (-z2 * theta - z3 * pw) at each time step
    z_tmp2 = zeros(N, 3 * N);
    for i = 1:N
        z_tmp2(i, (i-1) * 3 + 2) = -theta;
        z_tmp2(i, (i-1) * 3 + 3) = -pw;
    end
%tmp2 - coefficient to satisfy equality constraint: x - z2 * theta - z3 * P = 0
    tmp2 = [x_tmp, zeros(N), z_tmp2];
%z_tmp3 - sub-vector of Aeq to obtain (z1 + z2 + z3) at each time step
    z_tmp3 = zeros(N, 3 * N);
    for i = 1:N
        z_tmp3(i, (i-1)*3 + 1 : i * 3) = 1;
    end
%tmp3 - coefficient to satisfy equality constraint: z1 + z2 + z3 = 1
    tmp3 = [zeros(N, 2 * N), z_tmp3];
    Aeq = [Aeq; tmp1; tmp2; tmp3];
%beq - the vector for equality constraints
    beq = zeros(1);
    beq(1, 1) = e/t;
%beq_tmp1 - sub-vector of beq to satisfy equality constraints: y - z2 - z3 = 0 and x - z2 * small_const - z3 * P = 0
%beq_tmp2 - sub-vector of beq to satisfy equality constraints: z1 + z2 + z3 = 1
    beq_tmp1 = zeros(2 * N, 1);
    beq_tmp2 = ones(N, 1);
    beq = [beq;beq_tmp1;beq_tmp2];
    
%   Inequality Constraints - A*x <= b
%       A,b are the matrix and vector components respectively
%       Include upper and lower bounds for each variable

A_bound_up = eye(N);
A_bound_low = -1 * A_bound_up;

%A - the matrix for inequality constraints
%construct the sub-matrix of A corresponding to upper and lower bounds
    A = [A_bound_up; A_bound_low];
    A = [A,zeros(size(A))];
%b - the vector for inequality constraints
%b_bound_up - the vector correspond to upper bounds; b_bound_low - the vector correspond to lower bounds;
    b_bound_up = pw * ones(N, 1);
    b_bound_low = zeros(N, 1);
%construct the sub-vector of b corresponding to upper and lower bounds
    b = [b_bound_up; b_bound_low];

%   Ancillary minimum off-time constraints
%       construct the matrix for inequalities
    for i = 1: t_off - 1
%       construct the band matrix with elements 1; -1; (corresponding to the ON & OFF) and
%       1 on the i-th position (corresponding to switching-on after i intervals)
        first_col = [1;zeros(N - t_off - 1, 1)];
        first_row = [1, -1, zeros(1,i-1),1,zeros(1,N - 3 - i + 1)];

        band = toeplitz(first_col, first_row);
%supplement band with additional zeros corresponding to power status variables
        band = [zeros(size(band)), band];

        A = [A; band];
    end

    b_row = (N - t_off) * (t_off - 1);
    b = [b;ones(b_row, 1)]; %Linear inequality constraint vector
    
%   A_row - row number of A without considering z elements
    A_row = size(A, 1);
    A = [A, zeros(A_row, 3 * N)]; %Linear inequality constraint matrix

  lb = zeros(1,5*N); %Specifies the lower bound of all variables to be 0
  ub = [inf(1,N), ones(1,4*N)]; %Specifies no upper bound for x variables while y and z variables have an upper bound of 1
     
[solution,Final_Cost,~] = linprog(f',A,b,Aeq,beq,lb,ub);

%Specifies the portions of the solutions relating to x,y and z respectively
x = solution(1:N);
y = solution((N+1):(2*N));
z = solution((2*N)+1:5*N);

%Display final cost, power status and ON-OFF status
display(Final_Cost)
display(x);
display(y);

figure
%Plot of power status against time
subplot(2,1,1)
stairs(solution(1:N),'black')
xlabel('Time')
ylabel('Power Status (W)')
%Plot of price against time
subplot(2,1,2)
plot(pr,'black')
xlabel('Time')
ylabel('Price ($/Wh)')
end