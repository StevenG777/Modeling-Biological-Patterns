% MATH 132 Lecture Session

%% Discretize Space and Time
clc
close all
clear all

% The Time Step
delta_t = 0.1;
% The grid resolution
N = 100;
% The Step Size
delta_x = 1/(N-1);
% Help Sparse

% Let construct the matrix A
% Define it first as a sparse matrix
A_matrix = zeros(10,10);
A = sparse(A_matrix);

% Let fill the matrix
for i=1:N
    A_matrix(i,i) = -2 / delta_x^2;
    % If we are not on the left boundary add the left neighbor
    % contribution
    if(i>1)
         A_matrix(i,i-1) = 1 / delta_x^2;
    else
         A_matrix(i,N-1) = 1 / delta_x^2;
    end
    if(i<1)
         A_matrix(i,i+1) = 1 / delta_x^2;
    else
         A_matrix(i,2)   = 1 / delta_x^2;
    end
end

A_full = full(A_matrix);

%% Let's solve the heat equation
% initialize the solution
% this is the solution at the time step n
solution_n_vector =  zeros(N,1);
% this is the solution at the time step np1
solution_np1_vector =  zeros(N,1);

% maximum number of iterations
max_ite = 1000;

% let's initialize the solution vector
solution_n_vector = rand(N,1);

% Try first FWD
% We build the iteration matrix (I + delta * A)
% iteration_matrix = sparse(N,N);
iteration_matrix = A_matrix;
iteration_matrix = iteration_matrix * delta_t;
for i = 1:N
    iteration_matrix(i,i) = 1 + iteration_matrix(i,i);
end

for iteration = 1:max_ite
    % Compute the solution
    solution_np1_vector = iteration_matrix * solution_n_vector;
    figure(1)
    plot(solution_np1_vector)
    pause(0.01)
    % Update the solution: the new current step is the most recent one
    solution_n_vector = solution_np1_vector;
end

%{
% TPZ
Aplus = sparse(N,N);
Aplus = A_matrix;
Aplus = 0.5 * delta_t * Aplus;
for i = 1:N
    Aplus(i,i) = 1 + Aplus(i,i);
end


Aminus = sparse(N,N);
Aminus = A_matrix;
Aminus = -0.5 * delta_t * Aminus;
for i = 1:N
    Aminus(i,i) = 1 + Aminus;
end
%}
    