
clc
clear
A=delsq(numgrid('S',102));
L=ichol(A);
n=size(A,1);
b=A*ones(n,1);
tol=1.e-8;
maxit=50;


[x, flag, relres, iter, resvec] = pcg(A, b, tol, maxit, L, L');
display(flag); display(relres); display(iter);

[myx, myresvec, myiter] = mypcg(A, b, tol, maxit, L);
display(myiter);

semilogy(0:iter, resvec,'m-*', 0:myiter, myresvec,'b-o');
legend('PCG matlab','PCG implementation');
xlabel('Iteration ');
ylabel('Residual norm');