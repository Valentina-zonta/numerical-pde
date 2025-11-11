function [x, resvec, iter] = mycg(A, b, tol, maxit)

 % Input:
    % A: square symmetric positive definite matrix
    % b: right-hand side vector colum!!!
    % tol: tolerance for convergence
    % maxit: maximum number of iterations
    % L: Cholesky factorization of A (preconditioner)
   resvec = zeros(1,maxit);
    x(:,1)= zeros(1, size(A, 1));
    r(:,1) = b-A*x(:,1);
    p(:,1) = r(:,1);

   i = 1;


while (norm(r(:,i)) > tol * norm(b)) && (i < maxit)

    z(:,i) = A*p(:,i);
    alpha = (r(:,i)'*r(:,i))./(z(:,i)'*p(:,i));
    x(:,i+1) = x(:,i) + alpha*p(:,i);
    r(:,i+1) = r(:,i) - alpha*z(:,i);
    beta = (r(:,i+1)'*r(:,i+1))./(r(:,i)'*r(:,i));
    p(:,i+1) = r(:,i+1) + beta*p(:,i);
    resvec(i)=norm(r(:,i));
    i = i + 1;
end

resvec(i) = norm(b - A * x(:, i));
iter = i - 1;

