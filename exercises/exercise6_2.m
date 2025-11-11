

clc
clear
%set A
D = readmatrix('mat13041.rig.txt');
A = spconvert(D);
%set ILU
setup.type='crout';
setup.droptol=0.01;
[L,U]=ilu(A,setup);
%set B
i_values = 1:size(A,1);
i_values = i_values(:);
x_exact = 1 ./ sqrt(i_values);
b=A*x_exact;


%set imput elements
tol=10^(-10);
kmax=550;
x0=zeros(size(A,1),1);
restart=5000;

%solve with my precond. gmres
[x1, iter1, resvec1, flag1] = myprecgmres(A, b, tol, kmax, x0, L, U);
[x2,flag2, relres2,iter2, resvec2] = gmres(A, b,restart, tol, kmax, L, U, x0);

semilogy(0:iter2(2), resvec2,'m-*', 0:iter1, resvec1,'b-o');
legend('preconditioned gmres','myprecgmres');
xlabel('Iteration ');
ylabel('Residual norm');