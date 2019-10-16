function [ du_x ] = Dx_f( u )
%DX Linear operator approximating derivative of u w.r.t x
% Using finite difference forward approximation w/ periodic boundary
% conditions.
    
    [m,n]=size(u);
    h = (m-1)/(m-1);
    A = -speye(m*n);
    B = sparse(1:m*n-m,m+1:m*n,1,m*n,m*n);
    C = sparse(m*n-m+1:m*n,1:m,1,m*n,m*n);
    D=A+B+C;
    du_x= reshape((1/h)*D*u(:),[m,n]);


end