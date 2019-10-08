% --   FINITE ELEMENTS FOR THE 1D POISSON EQUATION   --

% - Caroline Myklebust - cmy002

close all;

n_vals = [10,20,40,80,160];
K = 2; 
u_true = @(x)x*(1-x);
tol = 10^-6;
xx = linspace(0,1,200);
h = zeros(length(n_vals),1);
Error_n = zeros(length(n_vals),1);

% Initializing plot
hold on;
plot(xx,u_truee(xx));
legend('True Solution');
xlabel('x');
ylabel('u(x)');

for i = 1:length(n_vals)
    
    n = n_vals(i);
    h(i) = 1/(n-1);
    x = linspace(0,1,n);
    x_inter = linspace(x(1)+(h(i)/2),x(n)-(h(i)/2),n-1);
    f = 2 * ones(n,1);
    
    %Setting up coefficient matrix
    e = ones(n,1);  
    A = (spdiags([-e 2*e -e], -1:1, n,n));
    
    % Boundary conditions, Setting Dirichlet on all boundaries
    A(1,1) = 1;
    A(1,2) = 0;
    A(n,n) = 1;
    A(n,n-1) = 0;
    f(1) = 0;
    f(n) = 0;
   
    % Using iterative solver on the linear system AU=b
    u_est = Steepest_descent((1/h(i))*A,f*h(i),tol);
    % Can alternatively solve directly using MATLABS Gaussian Elimination
    % method:
    %u_est = ((1/h(i)^2)*A)\f;
    
    Error_n(i) = 0;
    for j=1:n-1
        u_est_inter(j) = (u_est(j)+u_est(j+1))/2;
        e_new= (((u_true(x_inter(j))-u_est_inter(j))).^2);
        Error_n(i) = Error_n(i) + e_new;
    
    end
    Error_n(i)=sqrt(Error_n(i)/j);
 
    loglog(x,u_est);
    title(sprintf('True solution vs Estimated solution for n= %.0f',n_vals(i)));
    legend(sprintf('Estimated solution for n= %.0f',n_vals(i)));
    fprintf('The error for n = %.0f is %.6e \n',n_vals(i),Error_n(i));
    pause;
    
end

hold off;

% -- ERROR PLOT --

hold on;
figure;
new_plt = plot(n_vals,Error_n);
legend(new_plt,'Error(n)');
title('Error plot for estimated solution vs partition of interval');
xlabel('N partitions');
ylabel('Error');
hold off;
pause;

% -- ERROR TABLE PLOT --
disp('');
disp('         --  E R R O R   T A B L E   --');
disp('    _____________________________________');
disp('    h value    error      ratio     order');
disp('    _____________________________________');

for k=1:length(n_vals)
    if k==1
        error_table_=[h(k), Error_n(k), 0, 0];
    
    else
        ratio = Error_n(k-1)/Error_n(k);
        order = log(abs(ratio)) / log(abs(h(k-1)/h(k)));
        error_table_=[h(k), Error_n(k), ratio, order];
    end
    disp(error_table_)
end


function[utrue] = u_truee(x)

    if length(x) > 1
        
        utrue = zeros(length(x));
        for i = 1:length(x)
            utrue(i) = x(i)*(1-x(i));
        end
    else
        utrue = x*(1-x);
    end    
end

    