% Assignment 1, Q1 (a)
% 12/04/21

% LTE analysis.

close all
clear all
clc

% Set up variables.
syms x t u(x, t) h k x_hat t_hat gamma kappa theta tau

% Define the Taylor expansions.
u_taylor = taylor(u(x, t), [x, t], [x_hat, t_hat], 'Order', 5);

ujn = subs(u(x, t), [x, t], [x_hat, t_hat]);
ujp1np1 = subs(u_taylor, [x, t], [x_hat + h, t_hat + k]);
ujnp1 = subs(u_taylor, [x, t], [x_hat, t_hat + k]);
ujm1np1 = subs(u_taylor, [x, t], [x_hat - h, t_hat + k]);
ujm1n = subs(u_taylor, [x, t], [x_hat - h, t_hat]);
ujp1n = subs(u_taylor, [x, t], [x_hat + h, t_hat]);

% Defining the PDE:
% u_t - (kappa * u_xx) + (gamma * u) = 0

equation = gamma*u(x_hat, t_hat) + diff(u(x_hat, t_hat), t_hat) - kappa*diff(u(x_hat, t_hat), x_hat, x_hat);
equation_manipulated = (1/2) * diff(equation, t_hat) * k;

% Define the LTE.
lte = ((ujnp1 - ujn)/k) - ...
      (kappa*(ujm1n - (2*ujn) + ujp1n + ujm1np1 - (2*ujnp1) + ujp1np1)/(2*(h^2))) + ...
      (gamma*(((1 - theta)*ujn) + (theta*ujnp1)));


% Simplify and rearrange the LTE.
lte = lte - equation; % We know equation == 0, so LTE is unchanged.
lte = lte - equation_manipulated; % Ditto (this will become useful in a second. See note below)
lte = simplify(lte);
lte = collect(lte, k);
latex(lte)

% The output of the line above shows that if theta = 1/2, then the k term
% disappears and we have a method which is O(k^{2} + h^{2}), otherwise, the
% method is O(k + h^{2})

% NOTE: EXPLAINING SOME OF THE MANIPULATIONS.
% By playing around with our expression, we that k's only factor is the
% derivative of our equation WRT time divided by two, that is,
% 1/2 * d/dt(gamma*u - kappa*u_xx + u_tt), which, we know is also zero.
% Hence, we can subtract it away from our equation as well.
