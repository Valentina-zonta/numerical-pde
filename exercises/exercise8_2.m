clc
clear

% Load matrix from file
A = load('ML_laplace.mtx.txt');
A = spconvert(A);

% Plot settings
axes('YScale', 'log')
xlabel('Iterations ');
ylabel('Residual norm');
hold on

% GMRES parameters

droptolvalues = [2*10^(-2) 10^(-2) 3*10^(-3) 10^(-3) 10^(-4) 10^(-5)];
tol = 10^(-8);
maxit = 500;
x = ones(size(A,1), 1);
b = A * x;
restart = 50;



% Loop through droptol values
i = 1;
for mydroptol = droptolvalues

    tstart1 = tic;
    setup.type = 'crout';
    setup.droptol = mydroptol;
    [L, U] = ilu(A, setup);
    tprec = toc(tstart1);

    tstart2 = tic;
    [x, flag, relres, iter, resvec] = gmres(A, b, restart, tol, maxit, L, U);
    tsol = toc(tstart2);

    ttot = tprec + tsol;

    myiter = (iter(1) - 1) * restart + iter(2);

    ro = (nnz(L) + nnz(U) - 13041) / nnz(A);

    line(:, i) = [mydroptol, myiter, tprec, tsol, ttot, relres, ro];

    % Change legend names
    leg(i) = "Drop Tol: " + mydroptol;

    % Plot iteration vs. residual
    plot(0:myiter, resvec, '-o');

    i = i + 1;
end

% Display legend and plot
legend(leg(i - 6), leg(i - 5), leg(i - 4), leg(i - 3), leg(i - 2), leg(i - 1));
hold off

% Display results in a table
for i = 1:6
    table(:, i) = (line(:, i)');
end

display(table');
