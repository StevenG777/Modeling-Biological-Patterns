% MATH 132 Discussion Problem

%% Step 1: Discretize dx
N = 10;
a = 0;
b = 1;
dx = (b-a)/N;

b = zeros(N+1,1);

%% Step 2: Approximate u''(x)
u = zeros(N+1,1);
u(1) = 0;
u(2) = 0;

b = zeros(N+1,1);
hess_u = @(x,i) (u(i+1) - 2*u(i) + u(i+1))/dx^2;

%% Step 3: Define Au = b Form
A = zeros(N+1,N+1);
A(1,1) = 1;
A(end,end) = 1;
for i = 2:N
    A(i,i-1) = 1;
    A(i,i)   = -2;
    A(i,i+1) = 1;
end

%% Step 4: solve A\b = u
