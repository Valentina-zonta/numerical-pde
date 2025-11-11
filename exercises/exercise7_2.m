clear
clc
%set A
D = readmatrix('mat13041.rig.txt');
A = spconvert(D);
%set ILU
setup.type='crout';
setup.droptol=10^(-2);
[L,U]=ilu(A,setup);
%set b
n=size(A,1);
x_exact=1./sqrt(1:n);
b=A*x_exact';

%set imput elements
tol=10^(-12);
kmax=550;
x0=zeros(size(A,1),1);
restart=[10 20 30 50];

tStart = tic;
[x6, flag6, relres6, iter6, resvec6] = gmres(A, b,5000, tol, kmax, L, U, x0);
tEnd = toc(tStart);
CPUtime6 = tEnd;
iter6_tot=((iter6(1)-1)*5000)+iter6(2);



tStart = tic;
[x2, flag2, relres2, iter2, resvec2] = gmres(A, b, restart(2), tol, kmax, L, U, x0);
tEnd = toc(tStart);
CPUtime2 = tEnd;
iter2_tot=((iter2(1)-1)*restart(2))+iter2(2);

tStart = tic;
[x3, flag3, relres3, iter3, resvec3] = gmres(A, b, restart(3), tol, kmax, L, U, x0);
tEnd = toc(tStart);
CPUtime3 = tEnd;
iter3_tot=((iter3(1)-1)*restart(3))+iter3(2);

tStart = tic;
[x4, flag4, relres4, iter4, resvec4] = gmres(A, b, restart(4), tol, kmax, L, U, x0);
tEnd = toc(tStart);
CPUtime4 = tEnd;
iter4_tot=((iter4(1)-1)*restart(4))+iter4(2);


tStart = tic;
[x5, flag5, relres5, iter5, resvec5] = gmres(A, b,190, tol, kmax, L, U, x0);
tEnd = toc(tStart);
CPUtime5 = tEnd;
iter5_tot=((iter5(1)-1)*190)+iter5(2);

tStart = tic;
[x1, flag1, relres1, iter1, resvec1] = gmres(A, b, restart(1), tol, kmax, L, U, x0);
tEnd = toc(tStart);
CPUtime1 = tEnd;
iter1_tot=((iter1(1)-1)*restart(1))+iter1(2);

semilogy(0:iter1_tot, resvec1,'g-*',0:iter2_tot, resvec2,'m-*',0:iter3_tot, resvec3,'b-*',0:iter4_tot, resvec4,'r-*');
legend('restart 10','restart 20','restart 30','restart 40');
xlabel('Iteration ');
ylabel('Residual norm');


fprintf('iter');
display(iter1);
display(iter2);
display(iter3);
display(iter4);