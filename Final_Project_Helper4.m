% Guass Seidel with tolerance instead of max iters

format long;
clc;clear;close all;

% Setting up matrix system
A = [4 -1 -1 0 0 0; -1 4 -1 -1 0 0; 0 -1 4 -1 -1 0;...
     0 0 -1 4 -1 -1; 0 0 0 -1 4 -1; 0 0 0 0 -1 4];
b = [18;18;4;4;26;16];
x = [0;0;0;0;0;0];

% Tolerance 
maxerr = 1e-10;

% First error should be something huge
err = Inf;
gs_error = [];

% Initialize iteration #
iter = 0;

% Quick dimension check
[m,n] = size(A);
if m~=n
	error('Matrix A is not square')  
end

while(err > maxerr)
    x_old = x;
    
    for i = 1:n
        sum = 0;
        % first summation   
        for j = 1:i-1
            sum = sum + A(i,j)*x(j);   
        end
        % second summation + first one
        for j = i+1:n
            sum = sum + A(i,j)*x_old(j); 
        end
        % method update
        x(i) = (1/A(i,i)) * (b(i)-sum); 
    end
    % update iteration count
    iter = iter+1;
    sols(:,iter) = x;
    err = norm(A*sols(:,iter)-b);
    gs_error = [gs_error err];
end


fprintf('Method converged in %d iterations \n', iter);
disp(x);

figure;
semilogy(gs_error);