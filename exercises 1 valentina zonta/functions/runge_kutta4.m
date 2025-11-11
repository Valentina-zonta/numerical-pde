
function [outputVector, output2, lasty] = runge_kutta4(f, a, passo, maxt)
format long;
    h = passo;  % step size
    T = maxt;  % final time
    N = T / h;  % number of intervals
    

    % Initialize arrays
    t = 0:h:T;
    y = zeros(N+1, 1);
    y(1) = a;
    k1 = zeros(N+1, 1);
    k2 = zeros(N+1, 1);
    k3 = zeros(N+1, 1);
    k4 = zeros(N+1, 1);

    for i = 1:length(t)
        k1(i) = f(y(i));
        k2(i) = f(y(i) + h * k1(i) / 2);
        k3(i) = f(y(i) + h * k2(i) / 2);
        k4(i) = f(y(i) + h * k3(i));
        y(i+1) = y(i) + h * (k1(i) + 2 * k2(i) + 2 * k3(i) + k4(i)) / 6;
    end

    outputVector = y;
    output2 = y(2);
    lasty= y(N+1);
end
