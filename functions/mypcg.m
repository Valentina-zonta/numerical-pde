
function [x, resvec, iter] = mypcg(A, b, tol, maxit, L)

 % Input:
    % A: square symmetric positive definite matrix
    % b: right-hand side vector colum!!!
    % tol: tolerance for convergence
    % maxit: maximum number of iterations
    % L: Cholesky factorization of A (preconditioner)
   resvec = zeros(1,maxit);
    x(:,1)= zeros(1, size(A, 1));
    r(:,1) = b-A*x(:,1);

    v = L\r(:,1); 
    p(:,1) = L'\v;
    
    v = L\r(:,1);
    ro(:,1) = r(:,1)'*(L'\v);

   i = 1;
convergence = true;

while convergence && (i <= maxit)
    z(:, i) = A * p(:, i);
    alpha = ro(:, i) / (z(:, i)' * p(:, i));
    x(:, i + 1) = x(:, i) + alpha * p(:, i);
    r(:, i + 1) = r(:, i) - alpha * z(:, i);
    v = L \ r(:, i + 1);
    g(:, i + 1) = L' \ v;
    ro(:, i + 1) = r(:, i + 1)' * g(:, i + 1);
    beta = ro(:, i + 1) / ro(:, i);
    p(:, i + 1) = g(:, i + 1) + beta * p(:, i);
    resvec(i) = norm(r(:, i));
    
    convergence = resvec(i) > tol * norm(b);
    i = i + 1;
end

resvec(i) = norm(b - A * x(:, i));
iter = i - 1;

