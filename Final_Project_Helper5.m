%{
Problem 5e) Error Func
inputs: N, f, Utrue

* Create Laplacian
A  = Laplacian(N)
delta = 1/n-1
x = 0: delta: 1
y = 0: delta: 1 (y = x)

* Calculate RHS
F = zeros(N^2, 1)
for i = 1:N
    for j = 1:N
        k = i + (j - 1) * N
        F(k) =  f(x(i), y(j))
    end
end

Recall: idx = i + (j-1)*N

Utrue = zeros(N^2,1)
for i = 1:N
    for j = 1:N
        k = i + (j -1) * N
        Utrue (k) = Utrue (x(i), y(j))
    end
end

f = @(x,y) (sin(pi x ))^2 (sin(pi * y))^2 .....
Utrue = @(x,y)

error = norm(Utrue - F)
%}

%{
Problem 8
Need: u, v, f, g, dt, n
(u = sqrt(length(n))

Ru = zeros(N^2, 1)
Rv = " "

for i = 1:N
    for j = 1:N
        x = i + (j - 1) * N
        Ru(k) = u(k) + dt * f(u(k),v(k))
        Rv(k) = v(k) + dt * g(u(k),v(k))
    end
end

for z = 1:ts:tf
[Ru, Rv] = RHS()
end

u = M\Ru;
v = M\Rv;
%}