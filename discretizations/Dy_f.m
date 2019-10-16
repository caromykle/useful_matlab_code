function [ du_y ] = Dy_f(u)
% DY Linear operator approximating derivative of u 
% (which is assumed to be a matrix in 2 dimensions.) w.r.t y
% Using finite difference one sided approximation w/ periodic boundary
% conditions.

    [m,n]=size(u);
    h = (m-1)/(m-1);
    I=eye(m);
    e=ones(m,1);
    nn=zeros(m,1);
    T=spdiags([nn -1*e e],[-1 0 1],m,m);
    T(m,1) = 1;
    A = (kron(I,T));
    du_y = reshape((1/h)*A*u(:),[m,n]);


end