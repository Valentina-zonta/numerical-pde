function [y,CPUtime,num_steps,norm_error]=crank_nic_pcg(matrix_width,passo,matrix,T_final)
I=speye(matrix_width);
h=passo;
M=matrix;
N=T_final/h;

A=I+(h/2)*M;
B=I-(h/2)*M;
y_0 = ones(matrix_width,1);
y(:,1)=y_0;


tStart = tic;

for i=1:N

b=B*y(:,i);
y(:,i+1)=pcg(A,b,h^3,500);

end

tEnd = toc(tStart);


CPUtime = tEnd;
num_steps = N;


y_exact = load('accurate_solution');
error= abs(y(:,end)-y_exact);
norm_error=norm(error, Inf);
