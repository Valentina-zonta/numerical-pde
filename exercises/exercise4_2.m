clc
clear
tol=10^(-8);
maxit=10000;
A=gallery('wathen',100,100);
b=rand(30401,1);
M = sparse(diag(diag(A)));
L=ichol(A);
  
[x,flag,relres,iter,resvec] = pcg(A,b,tol, maxit); 
[x1, flag1, relres1, iter1, resvec1] = pcg(A, b, tol, maxit, M);
[x2, flag2, relres2, iter2, resvec2] = pcg(A, b, tol, maxit, L, L');

semilogy(0:iter, resvec,'m-o', 0:iter1, resvec1,'b-o',  0:iter2, resvec2,'c-o');
legend('no preconditoner','Jacobi','IC(0)');
xlabel('iteration number');
ylabel('Residual norm'); 


 