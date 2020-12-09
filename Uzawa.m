function [U, norm_res,norm_sol] = Uzawa(K,F, ndof_u, U0, omega, tol)
% written based on stationary iteration method sec 4.5.8
% of Numerical Methods in Matrix Computations (NMMC)
% M*U = B, where U = [y x]', B = [b c]', M = [H A;A' 0]
% |H  A|y   b
% |A' 0|x = c
% 
% define S = A'*inv(H)*A==> id S is SPD the following converges for
% 0 < omega < 2/lamda_max(S)
% H*y_(k+1) = b - A*x_k
%   x_(k+1) = x_k + omega*(A'*y_(k+1) -c)
%

H = K(1:ndof_u,1:ndof_u);           % K_uu
AT = K(ndof_u+1:end,1:ndof_u);      % K_pu
A = K(1:ndof_u,ndof_u+1:end);       % K_up

b = F(1:ndof_u);
c = F(ndof_u+1:end);

y = U0(1:ndof_u);
x = U0(ndof_u+1:end);

itr = 0;
norm_res = 1;
while norm_res > tol
    
   x_old = x;
   y_old = y;
   U_old = [y_old;x_old];
   
   y = H\(b-A*x_old);
   x = x_old + omega*(AT*y-c);
   
   U = [y;x];
   
   itr = itr+1;
   
   norm_res(itr) = norm(F-K*U);
   norm_sol(itr) = norm(U_old-U);
   
end

end