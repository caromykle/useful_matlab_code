% Solving a large sparse matrix by steepest descent
function[u_soln] = Steepest_descent(A,f,tol)

    [m,n] = size(A);
    
    % Setting the initial guess of u to be a zero vector
    u0 = zeros(1,n);

    err = 10;

    while err > tol
        res = f - A*u0;
        alpha = (transpose(res)*res)/(transpose(res)*A*res);
        u_new = u0 + alpha*res;
        
        %Re-initializing:
        u0 = u_new;
        err = res;
    end
    u_soln = u0;
end
