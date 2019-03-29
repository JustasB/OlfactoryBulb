function  v = RK4(f, h , y_in)
    %% Runge-Kutta implementation
    %% Inputs:
    % f function
    % h step
    % y_in initial value
    %% Output:
    % v solution
    
    theta = [0, 0.5, 0.5, 1];
    alfa = (1/6)*[1 2 2 1];

    k(1,1:2) = h*f(y_in);
    k(2,1:2) = h*f(y_in+theta(2).*k(1,:));
    k(3,1:2) = h*f(y_in+theta(3).*k(2,:));
    k(4,1:2) = h*f(y_in+k(3,:));
    
    v = y_in + [alfa*k(:,1), alfa*k(:,2)];
end