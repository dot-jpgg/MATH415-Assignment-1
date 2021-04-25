% Assignment 1, Q2.
% 17/04/21

close all
clear all
clc
set(0,'defaultTextInterpreter','latex'); % latex-ing

% Solving heat equation with reaction-diffusion.
t0 = 0;
tf = 5;

L = 10;

x0 = -L;
xf = L;

h = 1/50;
k = 1/50;

M = (xf - x0)/h;
N = (tf - t0)/k;

x = linspace(x0,xf,M+1);
t = linspace(t0,tf,N+1);

% Using PDE solver.
m = 0;
u = pdepe(m,@heatpde,@heatic,@heatbc,x,t);

% Plotting solution
figure(1)
imagesc(x,t,u)
colorbar
xlabel('$x$')
ylabel('$t$')
title(sprintf('$u_{t} = u_{xx} + 4u(1 - u^{2})$ on $%d \\le x \\le %d$ and $%d \\le t \\le %d$,',x0,xf,t0,tf))
print('using_pde_solver', '-dpng', '-r300');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solving using IMEX method.
p = k/(h^2);

% Specifying B.C.'s, I.C.'s and exact solution.
u_0 = @(x) heatic(x);

% Constructing LHS matrix
L_main = ones(1, M-1) * (1 + (2*p));
L_top = ones(1, M-2) * -p;
L_bottom = ones(1, M-2) * -p;

L = diag(L_bottom, -1) + diag(L_main, 0) + diag(L_top, 1);

% Adding initial condition
for j=1:M-1
    xj=(j*h) + x0;
    U(j,1) = u_0(xj);
end

% Solving system.
for n=1:N
    % Renaming for convenience.
    U_n=U(:,n);
    
    % Step (1)
    for j=1:M-1
        U_j=U_n(j);
        U_star(:,j)=U_j + (4*k*U_j*(1-(U_j^2)));
    end
    
    % Step (2)
    U(:,n+1)=sparse(L)\transpose(U_star);
end

% Creating new U matrix with B.C.'s and I.C.'s
U_new=zeros(M+1,N+1);
for j=1:M-1
    for n=1:N+1
        U_new(j+1,n)=U(j,n);
    end
end

% Adding boundary conditions
for n=2:N+1
    tn=(n-1)*k + t0;
    U_new(1,n)=g_0(tn);
    U_new(M+1,n)=g_1(tn);
end

% Adding initial condition
for j=2:M
    xj=((j-1)*h) + x0;
    U_new(j,1) = u_0(xj);
end

% Renaming and reformatting
U = transpose(U_new);

% Error analysis.
error_matrix=U - u;
error=norm(U(end,:)-u(end,:));

% Steady-state condition.
diff = norm(U(end-1,:)-U(end,:));
tol = 1e-6;
ss = diff < tol;

% Plotting approximation
figure(2)
imagesc(t,x,U)
colorbar
xlabel('$x$')
ylabel('$t$')
title(sprintf('$u_{t} = u_{xx} + 4u(1 - u^{2})$ on $%d \\le x \\le %d$ and $%d \\le t \\le %d$,',x0,xf,t0,tf))
print('numerical_approximation', '-dpng', '-r300');

% Plotting error
figure(3)
imagesc(t,x,error_matrix)
colorbar
xlabel('$x$')
ylabel('$t$')
title('Error in $u_{t} = u_{xx} + 4u(1 - u^{2})$')
print('error_in_numerical_approx', '-dpng', '-r300');


% Creating surface plot
figure(4)
s=surf(x,t,U);
s.EdgeColor = 'none';
xlabel('$x$')
ylabel('$t$')
zlabel('Temperature')
view(60,40)
title('Surface plot of $u_{t} = u_{xx} + 4u(1 - u^{2})$ solution')
print('surface_plot_of_approximation', '-dpng', '-r300');




% Defining functions for PDE solver.
function [c,f,s] = heatpde(x,t,u,dudx)
c = 1;
f = dudx;
s = (4*u) - (4*(u^3));
end

function u0 = heatic(x)
u0 = exp(-(x^2));
end

function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;
end