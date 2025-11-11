function [outputVectorfe, output2fe] = forward_Euler( a)
format long;
    h = 0.05;  % step size
    T = 6;  % final time
    N = T / h;  % number of intervals
    f=@(b) -5*b;
    
  % Initialize arrays
    t = 0:h:T;
    y = zeros(N+1, 1);
    y(1) = a;

    %implement  forward Euler
    for i=1:length(t)
        y(i+1)=y(i)+h*f(y(i));
    end    

    outputVectorfe=y;
    output2fe=y(2);