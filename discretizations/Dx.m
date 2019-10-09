function [ du_x ] = Dx( u )
% APPROXIMATING FIRST DERIVATIVE
%-------------------------------

% Constructing a linear operator matrix which returns a periodic, onesided
% approximation of the first partial derivative in x of a input matrix u.

    [m,n]=size(u);
    h = (m-1)/(m-1);
    A = speye(m*n);
    B = sparse(m+1:m*n, 1:m*n-m,-1,m*n,m*n);
    C = sparse(1:m,m*n-m+1:m*n,-1,m*n,m*n);
    D=A+B+C;
    du_x= reshape((1/h)*D*u(:),[m,n]);
    
    
end