function [x, iter, resvec, flag] = myprecgmres(A, b, tol, kmax, x0, L, U)


    r0 =(L*U)\( b - A * x0);
    k = 0;
    rho0 = norm(r0);
    beta = rho0;
    v(:,1) = r0 / beta;
    resvec(1)=rho0;
    flag=0;

    % GMRES iteration
    while  rho0 > tol * norm((L*U)\b) && k < kmax
        k = k + 1;
        z=A*v(:,k);
        w=L\z;
        v(:,k+1) = U\w;
         
        % Arnoldi process
       
        for j = 1:k
            h(j,k) = v(:,k+1)' * v(:,j);
            v(:,k+1) = v(:,k+1) - h(j,k) * v(:,j);
        end
         h(k+1,k) = norm(v(:,k+1));
        
        % Normalize vk1v
         v(:,k+1) =  v(:,k+1) /h(k+1,k);
        
        % Solve the least squares problem using QR factorization
       [Q,R] = qr(h);
      
        e1 = zeros(size(Q,1), 1);
        e1(1) = 1;

        g=beta*Q'* e1;

        %X = lsqminnorm(R,g)
        rho0=abs(beta*Q(1,k+1));
        display(rho0);
        resvec(k+1)=rho0;

        if rho0 == 0
            flag = -1;
        end
      
    end
    iter=k;
    y = lsqminnorm( R,g);
    x=x0+(U\(v(:,1:end-1)*y));



