function[A]=periodic_coeff(m)
% Constructing the FD coefficient matrix with periodic boundary conditions.
% Returns a matrix of size m*m x m*m
% This matrix A is used to solve for U using AU=b.
% NOTE: is NOT divided by h^2, do this manually after return.
I=eye(m);
e=ones(m,1);
T=spdiags([e -4*e e],[-1 0 1],m,m);
T(1,m)=1; T(m,1) = 1;
K=sparse(m,m);
K(1,m)=1;
K(m,1)=1;
S=spdiags([e e],[-1 1],m,m);
A = (kron(I,T)+kron(S,I)+kron(K,I));

end

