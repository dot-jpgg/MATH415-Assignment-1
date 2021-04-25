% Assignment 1, Q3.
% 18/04/21

close all
clear all
clc
set(0,'defaultTextInterpreter','latex'); % latex-ing

% Solving the Fitzhugh-Nagumo system.
R = -0.04; % R=0.04 for (a), R=0 for (b), R=-0.04 for (c).
rho = 0.3;
D = 0.04;

t0 = 0;
tf = 13;

L = 10;

x0 = -L;
xf = L;
y0 = -L;
yf = L;

h = 1/10; % dx = dy = h
k = 1/10;

M = (xf - x0)/h; % In this case, Mx = My = 20/h.
N = (tf - t0)/k;

x = linspace(x0,xf,M-1);
y = linspace(y0,yf,M-1);
t = linspace(t0,tf,N+1);

% Solving using method (4) in assignment.
qu=(k*D)/(h^2);
qv=k/(h^2);
Au=sparse(getA(M-1,qu));
Av=sparse(getA(M-1,qv));

% Specifying I.C.'s.
u_0 = @(x,y) u_ic(x,y);
v_0 = @(x,y) v_ic(x,y);

% Adding initial condition
for ii=1:M-1
    for jj=1:M-1
        yj=(jj*h) + y0;
        xi=(ii*h) + x0;
        U(jj,ii) = u_0(xi,yj); % U(y,x)=u_0(x,y)
        V(jj,ii) = v_0(xi,yj);
    end
end

% Solving system.
for n=1:N
    % Renaming for convenience.
    U_n=U(:,:,n); % n-th dimension of the matrix.
    V_n=V(:,:,n);

    % Step (1) ------------------------------------------------------------
    
    F = zeros(1, (M-1)^2);
    
    for ii=1:M-1
        for jj=1:M-1
            row=jj+((ii-1)*(M-1)); % What row are we in?
            
            % Boundary Conditions.
            if ii == 1
                F(row) = F(row) + qu;
            end
            if ii == M-1
                F(row) = F(row) + qu;
            end
            if jj == 1
                F(row) = F(row) + qu;
            end
            if jj == M-1
                F(row) = F(row) + qu;
            end
            
            U_ji=U_n(jj,ii);
            
            U_col(row,1)=U_n(jj,ii);
            V_col(row,1)=V_n(jj,ii);
            
            U_star(row,1)=U_ji + (k*R*(U_ji^2)) - (k*(U_ji^3));
            
        end
    end
    i=1:(M-1)^2;j=1:(M-1)^2;
    s=F;m=(M-1)^2;nn=(M-1)^2;
    %Fu = sparse(diag(F, 0));
    Fu = sparse(i,j,s,m,nn);
    
    % Step (2) ------------------------------------------------------------
    U_np1=(sparse(Au-Fu))\(U_star + (k*(U_star-R)) - (rho*k*(V_col-U_star)));
    
    F = zeros(1, (M-1)^2);
    
    for ii=1:M-1
        for jj=1:M-1
            row=jj+((ii-1)*(M-1));
            
            % Boundary Conditions.
            if ii == 1
                F(row) = F(row) + qv;
            end
            if ii == M-1
                F(row) = F(row) + qv;
            end
            if jj == 1
                F(row) = F(row) + qv;
            end
            if jj == M-1
                F(row) = F(row) + qv;
            end
            
            U_np1_mat(jj,ii)=U_np1(row,1);
        end
    end
    U(:,:,n+1)=U_np1_mat;
    %Fv = sparse(diag(F, 0));
    i=1:(M-1)^2;j=1:(M-1)^2;
    s=F;m=(M-1)^2;nn=(M-1)^2;
    Fv = sparse(i,j,s,m,nn);
    
    % Step (3) -----------------------------------------------------------
    V_np1=(sparse(Av-Fv))\(V_col - (k*(V_col-U_col)));
    
    for ii=1:M-1
        for jj=1:M-1
            row=jj+((ii-1)*(M-1));
            
            V_np1_mat(jj,ii)=V_np1(row,1);
        end
    end
    V(:,:,n+1)=V_np1_mat;
end

U0=U(:,:,1);
Ufm1=U(:,:,end-1);
Uf=U(:,:,end);
Ufy=U(:,end,end);

V0=V(:,:,1);
Vfm1=V(:,:,end-1);
Vf=V(:,:,end);

% Steady-state condition.
diff = max(norm(Uf-Ufm1),norm(Vf-Vfm1));
tol = 1e-6;
ss = diff < tol;

% Producing surface plots of steady state.
figure(1)
imagesc(x,y,Uf)
colorbar
title("Fitzhugh-Nagumo in Steady-State $(U_{f})$")
xlabel("$x$")
ylabel("$y$")
print('Uf', '-dpng', '-r300');

figure(2)
imagesc(x,y,Vf)
colorbar
title("Fitzhugh-Nagumo in Steady-State $(V_{f})$")
xlabel("$x$")
ylabel("$y$")
print('Vf', '-dpng', '-r300');

figure(3)
imagesc(x,y,U0)
colorbar
title("Fitzhugh-Nagumo Initial Condition $(U_{0})$")
xlabel("$x$")
ylabel("$y$")
print('U0', '-dpng', '-r300');

figure(4)
imagesc(x,y,V0)
colorbar
title("Fitzhugh-Nagumo Initial Condition $(V_{0})$")
xlabel("$x$")
ylabel("$y$")
print('V0', '-dpng', '-r300');

% Animate:
figure(5)
Un=U(:,:,1);
imagesc(x,y,Un)
colorbar
title(sprintf("Fitzhugh-Nagumo System ($U$ variable) ($t = %f$)", 0))
xlabel("$x$")
ylabel("$y$")
gif('U.gif', 'DelayTime', 0)
for n=2:N+1
    figure(5)
    Un=U(:,:,n);
    imagesc(x,y,Un)
    colorbar
    title(sprintf("Fitzhugh-Nagumo System ($U$ variable) ($t = %f$)", n*k))
    xlabel("$x$")
    ylabel("$y$")
    gif
end

figure(6)
Vn=V(:,:,1);
imagesc(x,y,Vn)
colorbar
title(sprintf("Fitzhugh-Nagumo System ($V$ variable) ($t = %f$)", 0))
xlabel("$x$")
ylabel("$y$")
gif('V.gif', 'DelayTime', 0)
for n=2:N+1
    figure(6)
    Vn=V(:,:,n);
    imagesc(x,y,Vn)
    colorbar
    title(sprintf("Fitzhugh-Nagumo System ($V$ variable) ($t = %f$)", n*k))
    xlabel("$x$")
    ylabel("$y$")
    gif
end

% Defining initial conditions.
function u0 = u_ic(x,y)
u0=1;
end

function v0 = v_ic(x,y)
v0=exp(-(x^2+y^2));
end

% Defining Laplacian FD matrix
function A = getA(M,q)
% Constructing T
T_main = ones(1, M) * (1+(4*q)); % for the case that dx = dy = h
T_top = ones(1, M - 1) * (-q);
T_bottom = ones(1, M - 1) * (-q);

T = sparse(diag(T_bottom, -1) + diag(T_main, 0) + diag(T_top, 1));

% Constructing I
I = (-q)*sparse(eye(M));

% Constructing A
A = cell(M);
for ii = 1:M
    for jj = 1:M
        A(ii,jj) = {zeros(M)};
    end
end

for p = 1:M
    if p ~= M % Creating bottom diag
        A(p + 1, p) = {sparse(I)};
    end
    A(p, p) = {sparse(T)}; % Creating main diag
    if p ~= 1
        A(p - 1, p) = {sparse(I)}; % Creating top diag
    end
end

A = sparse(cell2mat(A)); % Creating A
end
















% figure(1)
% imagesc(x,y,U0)
% colorbar
% xlabel('$x$')
% ylabel('$y$')
% 
% figure(2)
% imagesc(x,y,Uf)
% colorbar
% xlabel('$x$')
% ylabel('$y$')
% 
% figure(3)
% imagesc(x,y,V0)
% colorbar
% xlabel('$x$')
% ylabel('$y$')
% 
% figure(4)
% imagesc(x,y,Vf)
% colorbar
% xlabel('$x$')
% ylabel('$y$')
% 
% figure(5)
% t=linspace(t0,tf,length(Ufy));
% plot(t,Ufy')

