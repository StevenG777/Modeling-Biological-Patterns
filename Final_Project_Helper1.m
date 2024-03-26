% MATH 132 Discussion 
% Baixi Guo

n = 10;

u = randn(100,1);
figure(1);
imagesc(reshape(u,n,n));
title('u(x,y)');
xlabel('x')
ylabel('y')

figure(2);
v = randn(100,1);
imagesc(reshape(v,n,n));
title('v(x,y)');
xlabel('x');
ylabel('y');
