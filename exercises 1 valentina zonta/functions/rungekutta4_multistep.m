function rungekutta4_multistep(a,b,c,d,xzero,yzero,f,g , passo, maxt)



h = passo;  % step size
T = maxt;  % final time
N = T / h;  % number of intervals
t = 0:h:T;

x(1)=19;
y(1)=22;

f = @(x, y)  x*(a-(b*y));
g = @(x, y)  y*((c*x)-d);


for i=1:length(t)
     k1x(i) = f(x(i),y(i));
     k1y(i) = g(x(i),y(i));

     k2x(i) = f(x(i) + h * k1x(i) / 2, y(i) + h * k1y(i) / 2);
     k2y(i) = g(x(i) + h * k1x(i) / 2, y(i) + h * k1y(i) / 2);

     k3x(i) = f(x(i) + h * k2x(i) / 2, y(i) + h * k2y(i) / 2);
     k3y(i) = g(x(i) + h * k2x(i) / 2, y(i) + h * k2y(i) / 2);

     k4x(i) = f(x(i) + h * k3x(i), y(i) + h * k3y(i));
     k4y(i) = f(x(i) + h * k3x(i), y(i) + h * k3y(i));
        
     x(i+1) = x(i) + h * (k1x(i) + 2 * k2x(i) + 2 * k3x(i) + k4x(i)) / 6;
     y(i+1) = y(i) + h * (k1y(i) + 2 * k2y(i) + 2 * k3y(i) + k4y(i)) / 6;
end    