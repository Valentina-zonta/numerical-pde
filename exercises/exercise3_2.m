clc
clear



fin=(10^(4)-5);
% Define the diagonal values
diagonalValues = [200*(1:5), ones(1, fin)];

% Create a sparse matrix using spdiags
sparseMatrix = spdiags(diagonalValues', 0, 10^(4), 10^(4));


n = 10000;
tol=10^(-8);
b=rand(n,1);
maxit=5000;



[x,flag,relres,iter,resvec] = pcg(sparseMatrix, b, 10^(-8), 30);
display(flag);
realresvec=resvec;

[xmy, resvecmy, itermy] = mycg(sparseMatrix, b, tol, 30);
itermy;

resvecmy2= resvecmy(1:7);

semilogy(0:itermy, resvecmy2,'r-*');
legend('cg');
xlabel('iteration number');
ylabel('Residual norm');