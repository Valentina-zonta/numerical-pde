% Define the parameters
alpha = 0.2;
beta = 0.01;
gamma = 0.004;
delta = 0.07;

% Initial conditions and time settings
x0 = 19;
y0 = 22;

T = 300;
h = 1e-3;

% Define the differential equations
f = @(x, y) x * (alpha - beta * y);
g = @(x, y) y * (gamma * x - delta);

% Number of timesteps
N = T/ h;

% Initialize arrays to store results
x_vals = zeros(1, N + 1);
y_vals = zeros(1, N + 1);
x_vals(1) = x0;
y_vals(1) = y0;

% 4th order Runge-Kutta method
for n = 1:N
    K1x = f(x_vals(n), y_vals(n));
    K1y = g(x_vals(n), y_vals(n));

    K2x = f(x_vals(n) + h * K1x / 2, y_vals(n) + h * K1y / 2);
    K2y = g(x_vals(n) + h * K1x / 2, y_vals(n) + h * K1y / 2);

    K3x = f(x_vals(n) + h * K2x / 2, y_vals(n) + h * K2y / 2);
    K3y = g(x_vals(n) + h * K2x / 2, y_vals(n) + h * K2y / 2);

    K4x = f(x_vals(n) + h * K3x, y_vals(n) + h * K3y);
    K4y = g(x_vals(n) + h * K3x, y_vals(n) + h * K3y);

    x_vals(n + 1) = x_vals(n) + h * (K1x + 2 * K2x + 2 * K3x + K4x) / 6;
    y_vals(n + 1) = y_vals(n) + h * (K1y + 2 * K2y + 2 * K3y + K4y) / 6;
end

% Plot the numerical solution
t = 0:h:T;
figure;
plot(t, x_vals, 'b', t, y_vals, 'r');
xlabel('Time');
ylabel('Population');
title('Lotka-Volterra Predator-Prey Model');
legend('Prey (x)', 'Predator (y)');
