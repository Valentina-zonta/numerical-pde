function simspon2(c,d)
format long;
    h = 0.05;  % step size
    T = 6;  % final time
    N = T / h;  % number of intervals

    % Initialize arrays
    t = 0:h:T;
    y = zeros(N+1, 1);
    y(1) = 1;
    y(2)=c;
    
    

      % Perform the 2-step Simpson's method
     for i = 1:N-1
         y(i+2) = (y(i) + (h / 3) * (-5 * y(i) + 4 * (-5 * y(i+1))) )/(1+(5/3)*h );
     end

        % Plot the numerical solution
        %plot(t, y, 'b-', 'LineWidth', 2);
        %hold on;

        % create and plot the exact solution
        y_exact = zeros(size(t));
        for i = 1:length(t)
        y_exact(i) = exp(-5 * t(i));
        end

        %find the errors 
        error =zeros(size(t));
        for i=1:length(t)
            error(i)=y(i)-y_exact(i);
        end
        

        figure(d)
        plot(t,error, '-r', 'LineWidth', 2);
        hold on 

        %Set plot labels and legend
        xlabel('t');
        ylabel('error_y');
        title('error for dy/dt=-5y ,y(0) = 1, y(1)');
    
        hold off;

        %plot(t_exact, y_exact, 'r--', 'LineWidth', 2);

        % Set plot labels and legend
        %xlabel('t');
        %ylabel('y');
        %title('Solution of y'' = -5y, y(0) = 1');
        %legend('Numerical Solution', 'Exact Solution');
        %hold off;


