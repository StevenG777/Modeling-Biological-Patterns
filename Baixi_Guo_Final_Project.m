%% MATH 132 Final Project
% Baixi Guo
% 05/12/2022
% Coding Sections

%% Problem 4 (c) & (d)
% Part(c) 
% Implemented in Function Definition 1 & 2

% Part (d)
clc;
N = 5;
delta_x = 1/(N-1);
A = Laplacian(N,delta_x);
disp(A);

%% Problem 4 (e)
% Implemented in Function Definition 3 & 4
% Anonymous Function
clc;
RHS = @(x,y) (sin(pi*x))^2 * (sin(pi*y))^2 - 2*pi^2*(cos(2*pi*x) * (sin(pi*y))^2 + (sin(pi*x))^2 * cos(2*pi*y));
Solution = @(x,y) sin(pi*x)^2 * sin(pi*y)^2;

% All N values
N = [5, 10, 20]';
err = zeros(3,1);

% Exact, Numerical Solution & Error
for i = 1 : size(N)
    delta_x = 1/(N(i)-1);
    A = Laplacian(N(i),delta_x);
    U_exact = Exact_Solution(N(i),Solution);
    U = Linear_System(N(i),A,RHS);
    err(i) = norm(U-U_exact);
end
disp(err);

%% Problem 6 (c)
% Calculate Abs(Eigenvalue of M)
clc;
D = 10e-5;
N = 5;
delta_x = 1/(N-1);
delta_t = 0.1 * delta_x^2;
A = Laplacian(N, delta_x);
M = eye(N^2) - D * delta_t * A;
ei = (abs(eig(M)));
disp(ei);

% Calculate Abs(Eigenvalue of M^-1)
ei = (abs(eig(inv(M))));
disp(max(ei));

%% Problem 6 (d) (e)
clc;
% Part (d)
% Implementation in Function Definition 5

% part (e)
% First Set
N = 3;
D = 10e-3;
delta_x = 1/(N-1);
delta_t = 0.1 * delta_x^2;
M = New_M_Matrix(N,D,delta_t);
ei = (abs(eig(M)));
disp(ei);

% Second Set
N = 4;
D = 10e-4;
delta_x = 1/(N-1);
delta_t = 0.1 * delta_x^2;
M = New_M_Matrix(N,D,delta_t);
ei = (abs(eig(M)));
disp(ei);

% Third Set
N = 6;
D = 10e-6;
delta_x = 1/(N-1);
delta_t = 0.1 * delta_x^2;
M = New_M_Matrix(N,D,delta_t);
ei = (abs(eig(M)));
disp(ei);

%% Problem 8 
% Implementation in Function Definition 6

%% Problem 9 Part(a)
% Set up parameters
N = 100;
t0 = 0;
tf = 100;
delta_t = (tf - t0)/500;
F = 0.09;
k = 0.059;
D = 10e-5;

% Set Up Anonymous Functions
f = @(u,v) -u * v^2 + F * (1 - u);
g = @(u,v) u * v^2 - (F + k) * v;

% Uniform Initial Condition
U = 1 * ones(N^2,1);
V = 2 * ones(N^2,1);
PDE_System_Solver(U,V,t0,tf,delta_t,N,f,g,D);

% Random Initial Condition
U = rand(N^2,1);
V = 1 - U;
PDE_System_Solver(U,V,t0,tf,delta_t,N,f,g,D);

% Localized Initial Condition
n = randi([1 N^2]);
U = ones(N^2,1);
U(n) = 0;
V = 1 - U;
PDE_System_Solver(U,V,t0,tf,delta_t,N,f,g,D);

%% Problem 9 Part(b)
% Set up parameters
N = 100;
t0 = 0;
tf = 100;
delta_t = (tf - t0)/50;
F = 0.029;
k = 0.057;
D = 10e-5;

% Set Up Anonymous Functions
f = @(u,v) -u * v^2 + F * (1 - u);
g = @(u,v) u * v^2 - (F + k) * v;

% Random Initial Condition
U = rand(N^2,1);
V = 1 - U;
PDE_System_Solver(U,V,t0,tf,delta_t,N,f,g,D);

%% Problem 9 Part (c)
N = 100;
t0 = 0;
tf = 100;
delta_t = (tf - t0)/50;
F = 0.037;
k = 0.06;
D = 10e-5;

% Set Up Anonymous Functions
f = @(u,v) -u * v^2 + F * (1 - u);
g = @(u,v) u * v^2 - (F + k) * v;

% Random Initial Condition
U = rand(N^2,1);
V = 1 - U;
PDE_System_Solver(U,V,t0,tf,delta_t,N,f,g,D);


%% Function Definition
% Definition 1
function [n,nL,nR,nT,nB] = Final_Index(N,i,j)
    n = i + (j - 1) * N;
    
    if i > 1
        nL = i - 1 + (j - 1) * N;
    elseif i == 1
        nL = i + N - 2 + (j - 1) * N;
    end

    if i < N
        nR = i + 1 + (j - 1) * N;
    elseif i == N
        nR = i - N + 2 + (j - 1) * N;
    end

    if j < N
        nT = i + j * N;
    elseif j == N
        nT = i + (j - N + 1) * N;
    end

    if j > 1
        nB = i + (j - 2) * N;
    elseif j == 1
        nB = i + (j + N - 3) * N;
    end
end

% Definition 2
function A = Laplacian(N,delta_x)
    A = zeros(N,N);

    for i = 1:N
        for j = 1:N
            [n,nL,nR,nT,nB] = Final_Index(N,i,j);
            A(n,n) = -4/(delta_x^2);
            A(n,nL) = 1/(delta_x^2);
            A(n,nR) = 1/(delta_x^2);
            A(n,nT) = 1/(delta_x^2);
            A(n,nB) = 1/(delta_x^2);
        end
    end
end

% Definition 3
function U = Linear_System(N,A,RHS) 
    % Set up
    B = zeros(N^2,1);
    M = eye(N^2) - A;
    x = linspace(0,1,N);
    y = linspace(0,1,N);
    
    % Get the B Function
    for i = 1:N
        for j = 1:N
            n = i + (j - 1) * N;
            B(n) = RHS(x(i),y(j));
        end
    end

    % Get the U numerical solution
    U = M\B;
end

% Definition 4
function U_true = Exact_Solution(N,Solution)
    U_true = zeros(N^2,1);
    x = linspace(0,1,N);
    y = linspace(0,1,N);

    for i = 1:N
        for j = 1:N
            n = i + (j - 1) * N;
            U_true(n) = Solution(x(i),y(j));
        end
    end
end

% Definition 5
function M = New_M_Matrix(N,D,delta_t)
    delta_x = 1/(N-1);
    A = Laplacian(N, delta_x);
    M = eye(N^2) - D * delta_t * A;
end

% Definition 6
function [U,V] = PDE_System_Solver(U,V,t0,tf,delta_t,N,f,g,D)
    % Set Up M matrix and Ru & Rv
    M = New_M_Matrix(N,D,delta_t);
    Ru = zeros(N^2,1);
    Rv = zeros(N^2,1);
   
    % Solve For U(n+1)
    for t = t0:delta_t:tf
        for i = 1 : N
            for j = 1 : N
                n = i + (j - 1) * N;
                Ru(n) = U(n) + delta_t * f(U(n),V(n));
                Rv(n) = V(n) + delta_t * g(U(n),V(n));
            end
        end
        
        U = M\Ru;
        V = M\Rv;
        
        % Plotting
        figure(1);
        imagesc(reshape(U,N,N));
        title('U(x,y)')
        xlabel("x")
        ylabel("y")
        caxis([0 1])
        colorbar
    end
end
